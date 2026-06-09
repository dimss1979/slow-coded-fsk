#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include "scf_private.h"

#define SYM_PHASES 4
#define FFT_RATIO 4
#define TONE_SPAN 8
#define WEIGHT_SCALE (255.0f)
#define SYNC_RATIO (2.3f)
#define MAX_CFO 200 /* Hz */

#define FFT_LEN (SCF_BB_SYM_LEN * FFT_RATIO)
#define SYNC_THR (SYNC_RATIO * SCF_SYNC_LEN * WEIGHT_SCALE * (1.0f / (2.0f * TONE_SPAN)))
#define SYNC_MAX (SCF_SYNC_LEN * WEIGHT_SCALE * 1.0f)

struct sym_phase {
    complex float *source;
    uint8_t demod_buf[SCF_PKT_MAX][FFT_LEN];
};

static struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan;
static fftwf_complex *fft_buf;
static complex float input_signal[SCF_BB_SYM_LEN * 2];
static float carrier_freq;
static float carrier_phase;
static float fft_window[SCF_BB_SYM_LEN];
static size_t demod_buf_idx;
static unsigned int symbol_counter;
static unsigned int symbol_skip;

static const float cfo_step = (float) SCF_BB_SRATE / (float) FFT_LEN;
static const float cfo_low_tone = -((float) SCF_BB_SRATE / (float) SCF_BB_SYM_LEN) * ((float) (SCF_TONES - 1) / 2.0f);
static const float cfo_freq_min = cfo_low_tone - MAX_CFO;
static const float cfo_freq_max = cfo_low_tone + MAX_CFO;
static const int cfo_bin_min = (int) (cfo_freq_min / cfo_step) + FFT_LEN / 2;
static const int cfo_bin_max = (int) (cfo_freq_max / cfo_step) + FFT_LEN / 2;

static bool initialized;

static void downconvert(float *signal)
{
    memcpy(&input_signal[0], &input_signal[SCF_BB_SYM_LEN], SCF_BB_SYM_LEN * sizeof(input_signal[0]));

    complex float baseband[SCF_SYM_LEN];
    complex float baseband_filtered[SCF_SYM_LEN];
    complex float carrier = cexpf(-I * carrier_phase);
    complex float carrier_step = cexpf(-I * 2.0f * M_PI * carrier_freq * (1.0f / (float) SCF_SRATE));
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        baseband[i] = signal[i] * carrier;
        carrier *= carrier_step;
    }
    carrier_phase += 2.0f * M_PI * carrier_freq * ((float) SCF_SYM_LEN / (float) SCF_SRATE);
    while (carrier_phase > 2.0f * M_PI) {
        carrier_phase -= 2.0f * M_PI;
    }

    scf_filter_rx(baseband_filtered, baseband);
    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        input_signal[SCF_BB_SYM_LEN + i] = baseband_filtered[i * SCF_DEC_RATIO];
    }
}

static void demodulate(struct sym_phase *c)
{
    for (size_t i = 0; i < FFT_LEN; i++) {
        if (i < SCF_BB_SYM_LEN) {
            fft_buf[i] = c->source[i];
            fft_buf[i] *= fft_window[i];
        } else {
            fft_buf[i] = 0.0f;
        }
    }

    fftwf_execute_dft(fft_plan, fft_buf, fft_buf);

    float power[FFT_LEN];

    for (int i = 0; i < FFT_LEN; i++) {
        complex float bin = fft_buf[i];
        size_t shifted_idx = (i + FFT_LEN / 2) % FFT_LEN;
        power[shifted_idx] = crealf(bin) * crealf(bin) + cimagf(bin) * cimagf(bin);
    }

    for (int i = TONE_SPAN * FFT_RATIO; i < FFT_LEN - (TONE_SPAN * FFT_RATIO); i++) {
        float power_sum = 0.0f;

        for (int j = -TONE_SPAN; j < TONE_SPAN; j++) {
            int n = i + j * FFT_RATIO;
            power_sum += power[n];
        }

        c->demod_buf[demod_buf_idx][i] = WEIGHT_SCALE * power[i] / power_sum;
    }
}

static void find_sync(struct sym_phase *p, uint32_t *sync_vector, size_t packet_fec_len, float *sync_weight, size_t *sync_bin)
{
    uint32_t max_weight = 0;
    size_t max_bin = 0;

    size_t step = (SCF_FEC_LEN * packet_fec_len) / SCF_SYNC_LEN + 1;

    size_t pos[SCF_SYNC_LEN];
    for (size_t i = 0; i < SCF_SYNC_LEN; i++) {
        pos[i] = (demod_buf_idx + (i - SCF_SYNC_LEN + 1) * step + SCF_PKT_MAX) % SCF_PKT_MAX;
    }

    for (int b = cfo_bin_min; b <= cfo_bin_max; b += FFT_RATIO) {
        uint32_t weight = 0;
        for (size_t i = 0; i < SCF_SYNC_LEN; i++) {
            size_t t = b + sync_vector[i] * FFT_RATIO;
            weight += p->demod_buf[pos[i]][t];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b;
        }
    }

    int b0 = max_bin - FFT_RATIO * 2;
    int b1 = max_bin + FFT_RATIO * 2;
    for (int b = b0; b < b1; b++) {
        uint32_t weight = 0;
        for (size_t i = 0; i < SCF_SYNC_LEN; i++) {
            size_t t = b + sync_vector[i] * FFT_RATIO;
            weight += p->demod_buf[pos[i]][t];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b;
        }
    }

    if (max_weight > SYNC_THR) {
        *sync_weight = (float) max_weight / SYNC_MAX;
        *sync_bin = max_bin;
    }
}

static uint8_t ml_decode(size_t *positions, size_t phase, size_t bin, size_t codeword_size, uint32_t *code_table)
{
    uint32_t max_weight = 0;
    uint8_t max_symbol = 0;
    for (size_t i = 0; i < 256; i++) {
        uint32_t weight = 0;
        for (size_t s = 0; s < codeword_size; s++) {
            size_t b = bin + FFT_RATIO * code_table[i * codeword_size + s];
            weight += sym_phase[phase].demod_buf[positions[s]][b];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_symbol = i;
        }
    }

    return max_symbol;
}

static void fft_window_init(void)
{
    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        fft_window[i] = sinf((M_PI * i) / SCF_BB_SYM_LEN);
    }
}

void scf_rx_init(float freq)
{
    symbol_counter = 0;
    symbol_skip = 0;
    carrier_freq = freq;

    if (!initialized) {
        for (size_t i = 0; i < SYM_PHASES; i++) {
            sym_phase[i].source = &input_signal[SCF_BB_SYM_LEN - i * SCF_BB_SYM_LEN / SYM_PHASES];
        }

        fft_buf = fftwf_alloc_complex(FFT_LEN);
        assert(fft_buf);

        fft_plan = fftwf_plan_dft_1d(
            FFT_LEN,
            fft_buf,
            fft_buf,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        assert(fft_plan);

        fft_window_init();

        initialized = true;
    }
}

void scf_rx(scf_rx_result *result, float signal[SCF_SYM_LEN])
{
    memset(result, 0, sizeof(scf_rx_result));
    result->symbol_counter = symbol_counter;

    downconvert(signal);

    float max_sync_weight = 0;
    size_t max_sync_bin = 0;
    size_t max_sync_phase = 0;
    size_t packet_type = 0;

    for (size_t i = 0; i < SYM_PHASES; i++) {
        demodulate(&sym_phase[i]);

        for (size_t pt = 0; pt < SCF_PACKET_TYPES; pt++) {
            float sync_weight = 0;
            size_t sync_bin = 0;
            find_sync(
                &sym_phase[i],
                &scf_packet_sync_vector[pt][0],
                scf_packet_fec_len[pt],
                &sync_weight,
                &sync_bin
            );

            if (sync_weight > max_sync_weight) {
                max_sync_weight = sync_weight;
                max_sync_bin = sync_bin;
                max_sync_phase = i;
                packet_type = pt;
            }
        }
    }

    if (max_sync_weight > 0) {

        int sync_bin_shifted = max_sync_bin - FFT_LEN / 2;
        float sync_cfo = sync_bin_shifted * cfo_step - cfo_low_tone;

        result->got_sync = true;
        result->sync_phase = max_sync_phase;
        result->sync_packet_type = packet_type;
        result->sync_cfo = sync_cfo;
        result->sync_weight = max_sync_weight;

        if (!symbol_skip) {
            size_t outer_codeword_len = scf_packet_fec_len[packet_type];

            size_t positions[SCF_MSG_FEC_MAX][SCF_FEC_LEN];

            size_t block_len = (SCF_FEC_LEN * outer_codeword_len) / SCF_SYNC_LEN;
            size_t first_block_len = block_len + (SCF_FEC_LEN * outer_codeword_len) % SCF_SYNC_LEN;
            size_t pos = (demod_buf_idx - (SCF_FEC_LEN * outer_codeword_len + SCF_SYNC_LEN) + 1 + SCF_PKT_MAX) % SCF_PKT_MAX;
            size_t sync_cnt = first_block_len;
            size_t data_cnt = 0;

            while (pos != demod_buf_idx) {
                if (sync_cnt) {
                    size_t inner_i = data_cnt / outer_codeword_len;
                    size_t outer_i = data_cnt % outer_codeword_len;

                    positions[outer_i][inner_i] = pos;
                    data_cnt++;
                    sync_cnt--;
                } else {
                    sync_cnt = block_len;
                }

                pos = (pos + 1) % SCF_PKT_MAX;
            }

            uint8_t outer_codeword[SCF_MSG_FEC_MAX];

            for (size_t i = 0; i < outer_codeword_len; i++) {
                uint8_t byte_val = ml_decode(
                    &positions[i][0],
                    max_sync_phase,
                    max_sync_bin,
                    SCF_FEC_LEN,
                    &scf_inner_code[0][0]
                );
                outer_codeword[i] = byte_val;
            }

            if (scf_packet_decode(result, outer_codeword)) {
                symbol_skip = 5;
            }
        }
    }

    demod_buf_idx = (demod_buf_idx + 1) % SCF_PKT_MAX;

    symbol_counter++;
    if (symbol_skip) {
        symbol_skip--;
    }
}
