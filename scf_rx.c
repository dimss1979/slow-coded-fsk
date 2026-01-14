#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <fec.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_packet.h"
#include "scf_rx.h"

#define SYM_PHASES 4
#define FFT_RATIO 4
#define CW_FILTER_LEN 10
#define DEC_RATIO 4
#define TONE_SPAN 8

#define BB_SRATE (SCF_SRATE / DEC_RATIO)
#define BB_SYM_LEN (SCF_SYM_LEN / DEC_RATIO)
#define FFT_LEN (BB_SYM_LEN * FFT_RATIO)

struct sym_phase {
    complex float *source;
    uint8_t demod_buf[FFT_LEN][SCF_PKT_MAX];
    uint8_t peak_buf[FFT_LEN][SCF_PKT_MAX];
};

typedef enum {
    S_PREAMBLE,
    S_DATA,
} e_state;

static struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan;
static fftwf_complex *fft_buf;
static complex float input_signal[BB_SYM_LEN * 2];
static float carrier_freq;
static float carrier_phase;
static complex float fir_tail[SCF_FIR_LEN_RF];
static float fft_window[BB_SYM_LEN];
static size_t demod_buf_idx;
static unsigned int symbol_counter;
static e_state state = S_PREAMBLE;
static size_t rx_bin;
static size_t rx_phase;
static uint32_t last_preamble_weight;
static size_t data_raw_len;
static size_t data_fec_len;
static size_t data_count;
static uint8_t data_fec_buf[SCF_MSG_FEC_MAX];
static uint32_t data_weight_buf[SCF_MSG_FEC_MAX];
static void *data_rs_code;

static bool initialized;

static void downconvert(float *signal)
{
    memcpy(&input_signal[0], &input_signal[BB_SYM_LEN], BB_SYM_LEN * sizeof(input_signal[0]));

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
    for (size_t i = 0; i < BB_SYM_LEN; i++) {
        input_signal[BB_SYM_LEN + i] = baseband_filtered[i * DEC_RATIO];
    }
}

static void demodulate(struct sym_phase *c)
{
    for (size_t i = 0; i < FFT_LEN; i++) {
        if (i < BB_SYM_LEN) {
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
        power[i] = crealf(bin) * crealf(bin) + cimagf(bin) * cimagf(bin);
    }

    for (int i = 0; i < FFT_LEN; i++) {
        float power_sum = 0.0f;
        uint8_t is_local_peak = 1;

        for (int j = -TONE_SPAN; j < TONE_SPAN; j++) {
            int n = (i + j * FFT_RATIO + FFT_LEN) % FFT_LEN;
            power_sum += power[n];

            if (power[n] > power[i]) {
                is_local_peak = 0;
            }
        }

        c->demod_buf[i][demod_buf_idx] = 255.0f * power[i] / power_sum;
        c->peak_buf[i][demod_buf_idx] = is_local_peak;
    }
}

static void find_preamble(struct sym_phase *p, uint32_t *sync_vector, uint32_t *preamble_weight, uint32_t *preamble_symbols, size_t *preamble_bin)
{
    uint32_t max_weight = 0;
    size_t max_bin = 0;

    for (int b = 0; b < FFT_LEN; b += FFT_RATIO) {
        uint32_t weight = 0;
        for (size_t i = 0; i < SCF_PREAMBLE; i++) {
            size_t t = (b + sync_vector[i] * FFT_RATIO) % FFT_LEN;
            size_t pos = (demod_buf_idx + i - SCF_PREAMBLE + SCF_PKT_MAX + 1) % SCF_PKT_MAX;
            weight += p->demod_buf[t][pos];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b;
        }
    }

    size_t b0 = max_bin - FFT_RATIO * 2;
    size_t b1 = max_bin + FFT_RATIO * 2;
    for (size_t b = b0; b < b1; b++) {
        size_t b_wrapped = (b + FFT_LEN) % FFT_LEN;
        uint32_t weight = 0;
        for (size_t i = 0; i < SCF_PREAMBLE; i++) {
            size_t t = (b + sync_vector[i] * FFT_RATIO) % FFT_LEN;
            size_t pos = (demod_buf_idx + i - SCF_PREAMBLE + SCF_PKT_MAX + 1) % SCF_PKT_MAX;
            weight += p->demod_buf[t][pos];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b_wrapped;
        }
    }

    uint32_t decoded_count = 0;
    for (size_t i = 0; i < SCF_PREAMBLE; i++) {
        size_t t = (max_bin + sync_vector[i] * FFT_RATIO) % FFT_LEN;
        size_t pos = (demod_buf_idx + i - SCF_PREAMBLE + SCF_PKT_MAX + 1) % SCF_PKT_MAX;

        if (p->peak_buf[t][pos]) {
            decoded_count++;
        }
    }

    if (decoded_count > SCF_PREAMBLE_THR) {
        *preamble_weight = max_weight;
        *preamble_symbols = decoded_count;
        *preamble_bin = max_bin;
    }
}

static uint8_t ml_decode(uint32_t *weight, size_t *positions, size_t codeword_size, uint32_t *code_table)
{
    uint32_t max_weight = 0;
    uint8_t max_symbol = 0;
    for (size_t i = 0; i < 256; i++) {
        uint32_t weight = 0;
        for (size_t s = 0; s < codeword_size; s++) {
            size_t b = (rx_bin + FFT_RATIO * code_table[i * codeword_size + s]) % FFT_LEN;
            weight += sym_phase[rx_phase].demod_buf[b][positions[s]];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_symbol = i;
        }
    }

    *weight = max_weight;
    return max_symbol;
}

static bool msg_decode(uint8_t *msg)
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
    for (size_t i = 0; i < BB_SYM_LEN; i++) {
        fft_window[i] = sin(M_PI * i / BB_SYM_LEN);
    }
}

void scf_rx_init(float freq)
{
    state = S_PREAMBLE;
    symbol_counter = 0;
    carrier_freq = freq;

    if (!initialized) {
        for (size_t i = 0; i < SYM_PHASES; i++) {
            sym_phase[i].source = &input_signal[BB_SYM_LEN - i * BB_SYM_LEN / SYM_PHASES];
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

size_t scf_rx_symbol(uint8_t *msg, float *signal)
{
    size_t msg_len = 0;

    downconvert(signal);

    uint32_t max_preamble_weight = 0;
    uint32_t max_preamble_symbols = 0;
    size_t max_preamble_bin = 0;
    size_t max_preamble_phase = 0;
    size_t packet_type = 0;

    for (size_t i = 0; i < SYM_PHASES; i++) {
        demodulate(&sym_phase[i]);

        for (size_t pt = 0; pt < SCF_PACKET_TYPES; pt++) {
            uint32_t preamble_weight = 0;
            uint32_t preamble_symbols = 0;
            size_t preamble_bin = 0;
            find_preamble(
                &sym_phase[i],
                &scf_packet_sync_vector[pt][0],
                &preamble_weight,
                &preamble_symbols,
                &preamble_bin
            );

            if (preamble_weight > max_preamble_weight) {
                max_preamble_weight = preamble_weight;
                max_preamble_symbols = preamble_symbols;
                max_preamble_bin = preamble_bin;
                max_preamble_phase = i;
                packet_type = pt;
            }
        }
    }

    bool preamble_found = false;
    if (max_preamble_weight > 0) {
        preamble_found = true;
    }

    switch (state) {
        case S_PREAMBLE:
            if (preamble_found) {
                printf(
                    " +++ Preamble at %i weight %u bin %li phase %li match %u/%u, packet of %li bytes\n",
                       symbol_counter, max_preamble_weight, max_preamble_bin, max_preamble_phase, max_preamble_symbols,
                       SCF_PREAMBLE, scf_packet_raw_len[packet_type]
                );
                last_preamble_weight = max_preamble_weight;
                rx_phase = max_preamble_phase;
                rx_bin = max_preamble_bin;
                data_raw_len = scf_packet_raw_len[packet_type];
                data_fec_len = scf_packet_fec_len[packet_type];
                data_rs_code = scf_packet_rs_code[packet_type];
                data_count = 0;
                state = S_DATA;
            }

            break;

        case S_DATA:
            if (data_count == 0 && preamble_found && max_preamble_weight > last_preamble_weight) {
                printf(
                    " +++ Stronger preamble at %i weight %u bin %li phase %li match %u/%u, packet of %li bytes\n",
                        symbol_counter, max_preamble_weight, max_preamble_bin, max_preamble_phase, max_preamble_symbols,
                       SCF_PREAMBLE, scf_packet_raw_len[packet_type]
                );
                last_preamble_weight = max_preamble_weight;
                rx_phase = max_preamble_phase;
                rx_bin = max_preamble_bin;
                data_raw_len = scf_packet_raw_len[packet_type];
                data_fec_len = scf_packet_fec_len[packet_type];
                data_rs_code = scf_packet_rs_code[packet_type];
                data_count = 0;

                break;
            }

            data_count++;

            if (data_count == SCF_FEC_LEN * data_fec_len) {
                for (size_t i = 0; i < data_fec_len; i++) {
                    size_t positions[SCF_FEC_LEN];
                    for (size_t j = 0; j < SCF_FEC_LEN; j++) {
                        size_t pos = (demod_buf_idx - data_count + 1 + i + j * data_fec_len + SCF_PKT_MAX) % SCF_PKT_MAX;
                        positions[j] = pos;
                    }

                    uint32_t byte_weight;
                    uint8_t byte_val = ml_decode(&byte_weight, positions, SCF_FEC_LEN, &scf_data_code[0][0]);
                    data_fec_buf[i] = byte_val;
                    data_weight_buf[i] = byte_weight;
                }

                if (msg_decode(msg)) {
                    msg_len = data_raw_len;
                }
                printf(" +++ Data finished at %i\n", symbol_counter);
                state = S_PREAMBLE;
            }

            break;
    }


    demod_buf_idx = (demod_buf_idx + 1) % SCF_PKT_MAX;

    symbol_counter++;

    return msg_len;
}
