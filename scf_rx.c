#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <fec.h>
#include <stdio.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_packet.h"
#include "scf_rx.h"

#define SYM_PHASES 8
#define FFT_RATIO 4

#define FFT_LEN (SCF_BB_SYM_LEN * FFT_RATIO)

struct sym_phase {
    complex float *source;
    float prev_peak;
    size_t prev_peak_idx;
    float tone_buf[SCF_TONES][SCF_PKT_MAX];
};

static struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan, ifft_plan;
static fftwf_complex *fft_buf;
static complex float input_signal[SCF_BB_SYM_LEN * 2];
static float carrier_freq;
static float carrier_phase;
static complex float fir_tail[SCF_FIR_LEN_RF];
static float fft_window[SCF_BB_SYM_LEN];
static float spectrum_sample[FFT_LEN];
static size_t demod_buf_idx;
static unsigned int symbol_counter;
static unsigned int symbol_skip;

static bool initialized;

static void downconvert(float *signal)
{
    memcpy(&input_signal[0], &input_signal[SCF_BB_SYM_LEN], SCF_BB_SYM_LEN * sizeof(input_signal[0]));

    complex float baseband[SCF_SYM_LEN];
    complex float baseband_filtered[SCF_SYM_LEN];
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        float carrier_i = sinf(carrier_phase);
        float carrier_q = cosf(carrier_phase);
        baseband[i] = signal[i] * carrier_i + signal[i] * I * carrier_q;

        carrier_phase += 2.0f * M_PI * carrier_freq * (1.0f / (float) SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }
    }
    scf_filter_rf(baseband_filtered, baseband, fir_tail);
    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        input_signal[SCF_BB_SYM_LEN + i] = baseband_filtered[i * SCF_DEC_RATIO];
    }
}

static void calculate_spectrum(float *spectrum, complex float *signal)
{
    for (size_t i = 0; i < FFT_LEN; i++) {
        if (i < SCF_BB_SYM_LEN) {
            fft_buf[i] = signal[i];
            fft_buf[i] *= fft_window[i];
        } else {
            fft_buf[i] = 0.0f;
        }
    }

    fftwf_execute(fft_plan);

    float sum = 0.0f;
    for (size_t i = 0; i < FFT_LEN; i++) {
        float magnitude = sqrtf(fft_buf[i] * conjf(fft_buf[i]));
        sum += magnitude;
        spectrum[i] = magnitude;
    }

    float dc_offset = sum / FFT_LEN;
    float energy = 1e-10f;
    for (size_t i = 0; i < FFT_LEN; i++) {
        spectrum[i] -= dc_offset;
        energy += spectrum[i] * spectrum[i];
    }

    float scale = sqrtf(energy);
    for (size_t i = 0; i < FFT_LEN; i++) {
        spectrum[i] /= scale;
    }
}

static void calculate_correlation(
    float *correlation, float *signal1, float *signal2
)
{
    complex float spectrum1[FFT_LEN];
    complex float spectrum2[FFT_LEN];

    for (size_t i = 0; i < FFT_LEN; i++)
        fft_buf[i] = signal1[i];

    fftwf_execute(fft_plan);

    for (size_t i = 0; i < FFT_LEN; i++)
        spectrum1[i] = fft_buf[i];

    for (size_t i = 0; i < FFT_LEN; i++)
        fft_buf[i] = signal2[i];

    fftwf_execute(fft_plan);

    for (size_t i = 0; i < FFT_LEN; i++)
        spectrum2[i] = fft_buf[i];

    for (size_t i = 0; i < FFT_LEN; i++)
        fft_buf[i] = spectrum1[i] * conjf(spectrum2[i]);

    fftwf_execute(ifft_plan);

    for (size_t i = 0; i < FFT_LEN; i++)
        correlation[i] = crealf(fft_buf[i]) / FFT_LEN;
}

static void demodulate(struct sym_phase *c)
{
    float received_spectrum[FFT_LEN];
    calculate_spectrum(received_spectrum, c->source);

    float correlation[FFT_LEN];
    calculate_correlation(correlation, spectrum_sample, received_spectrum);

    for (size_t t = 0; t < SCF_TONES; t++) {
        size_t tone_offset = (t + 1) * FFT_RATIO;
        size_t new_peak_idx = (c->prev_peak_idx + tone_offset) % FFT_LEN;
        float tone_weight = correlation[new_peak_idx] * c->prev_peak;
        c->tone_buf[t][demod_buf_idx] = tone_weight;
    }

    float peak = 0.0f;
    size_t peak_idx = 0;
    for (size_t i = 0; i < FFT_LEN; i++) {
        if (correlation[i] > peak) {
            peak = correlation[i];
            peak_idx = i;
        }
    }

    c->prev_peak = peak;
    c->prev_peak_idx = peak_idx;
}

static void find_sync(struct sym_phase *p, uint32_t *sync_vector, size_t packet_fec_len, float *sync_weight)
{
    size_t step = (SCF_FEC_LEN * packet_fec_len) / SCF_SYNC_LEN + 1;

    float weight = 0.0f;
    for (size_t i = 0; i < SCF_SYNC_LEN; i++) {
        size_t t = sync_vector[i];
        size_t pos = (demod_buf_idx + (i - SCF_SYNC_LEN + 1) * step + SCF_PKT_MAX) % SCF_PKT_MAX;
        weight += p->tone_buf[t][pos];
    }
    weight /= SCF_SYNC_LEN;

    if (weight > 0.01f)
        *sync_weight = weight;
}

static uint8_t ml_decode(size_t *positions, size_t phase, size_t codeword_size, uint32_t *code_table)
{
    float max_weight = 0.0f;
    uint8_t max_symbol = 0;
    for (size_t i = 0; i < 256; i++) {
        float weight = 0.0f;
        for (size_t s = 0; s < codeword_size; s++) {
            size_t tone_idx = code_table[i * codeword_size + s];
            weight += sym_phase[phase].tone_buf[tone_idx][positions[s]];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_symbol = i;
        }
    }

    return max_symbol;
}

static bool msg_decode(uint8_t *msg, uint8_t *data_fec_buf, void *data_rs_code, size_t data_raw_len)
{
    int symbol_error_count = decode_rs_char(data_rs_code, data_fec_buf, NULL, 0);

    if (symbol_error_count >= 0) {
        for (size_t i = 0; i < data_raw_len; i++) {
            data_fec_buf[i] ^= scf_data_scrambler[i];
        }
        printf(" +++ Outer FEC errors: %i\n", symbol_error_count);
        memcpy(msg, data_fec_buf, data_raw_len);
        return true;
    }

    return false;
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

        ifft_plan = fftwf_plan_dft_1d(
            FFT_LEN,
            fft_buf,
            fft_buf,
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );
        assert(ifft_plan);

        fft_window_init();

        initialized = true;
    }

    calculate_spectrum(spectrum_sample, &scf_waveform[0][0]);
}

size_t scf_rx_symbol(uint8_t *msg, float *signal)
{
    size_t msg_len = 0;

    downconvert(signal);

    float max_sync_weight = 0.0f;
    size_t max_sync_phase = 0;
    size_t packet_type = 0;

    for (size_t i = 0; i < SYM_PHASES; i++) {
        demodulate(&sym_phase[i]);

        for (size_t pt = 0; pt < SCF_PACKET_TYPES; pt++) {
            float sync_weight = 0.0f;
            find_sync(
                &sym_phase[i],
                &scf_packet_sync_vector[pt][0],
                scf_packet_fec_len[pt],
                &sync_weight
            );

            if (sync_weight > max_sync_weight) {
                max_sync_weight = sync_weight;
                max_sync_phase = i;
                packet_type = pt;
            }
        }
    }

    if (max_sync_weight > 0.0f) {
        printf(
            " +++ sync at %i weight %f phase %li, packet of %li bytes\n",
                symbol_counter, max_sync_weight, max_sync_phase,
                scf_packet_raw_len[packet_type]
        );

        if (!symbol_skip) {
            size_t data_raw_len = scf_packet_raw_len[packet_type];
            size_t data_fec_len = scf_packet_fec_len[packet_type];
            void *data_rs_code = scf_packet_rs_code[packet_type];

            size_t positions[SCF_MSG_FEC_MAX][SCF_FEC_LEN];

            size_t block_len = (SCF_FEC_LEN * data_fec_len) / SCF_SYNC_LEN;
            size_t first_block_len = block_len + (SCF_FEC_LEN * data_fec_len) % SCF_SYNC_LEN;
            size_t pos = (demod_buf_idx - (SCF_FEC_LEN * data_fec_len + SCF_SYNC_LEN) + 1 + SCF_PKT_MAX) % SCF_PKT_MAX;
            size_t sync_cnt = first_block_len;
            size_t data_cnt = 0;

            while (pos != demod_buf_idx) {
                if (sync_cnt) {
                    size_t inner_i = data_cnt / data_fec_len;
                    size_t outer_i = data_cnt % data_fec_len;

                    positions[outer_i][inner_i] = pos;
                    data_cnt++;
                    sync_cnt--;
                } else {
                    sync_cnt = block_len;
                }

                pos = (pos + 1) % SCF_PKT_MAX;
            }

            uint8_t data_fec_buf[SCF_MSG_FEC_MAX];

            for (size_t i = 0; i < data_fec_len; i++) {
                uint8_t byte_val = ml_decode(
                    &positions[i][0],
                    max_sync_phase,
                    SCF_FEC_LEN,
                    &scf_data_code[0][0]
                );
                data_fec_buf[i] = byte_val;
            }

            if (msg_decode(msg, data_fec_buf, data_rs_code, data_raw_len)) {
                msg_len = data_raw_len;
                symbol_skip = 5;
            }
        }
    }

    demod_buf_idx = (demod_buf_idx + 1) % SCF_PKT_MAX;

    symbol_counter++;
    if (symbol_skip) {
        symbol_skip--;
    }

    return msg_len;
}
