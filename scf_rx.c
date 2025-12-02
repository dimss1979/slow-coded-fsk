#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_packet.h"
#include "scf_rx.h"

#define SYM_PHASES 4
#define FFT_RATIO 4
#define CW_FILTER_LEN 13
#define DEC_RATIO 4

#define BB_SRATE (SCF_SRATE / DEC_RATIO)
#define BB_SYM_LEN (SCF_SYM_LEN / DEC_RATIO)
#define FFT_LEN (BB_SYM_LEN * FFT_RATIO)

struct sym_phase {
    complex float *source;
    float cw_filter_buf[FFT_LEN][CW_FILTER_LEN];
    int8_t demod_buf[FFT_LEN][SCF_PREAMBLE][SCF_TONES];
};

typedef enum {
    S_PREAMBLE,
    S_HEADER,
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
static uint8_t preamble_ref[SCF_PREAMBLE];
static size_t demod_buf_idx;
static unsigned int symbol_counter;
static e_state state = S_PREAMBLE;
static size_t rx_bin;
static size_t rx_phase;
static size_t header_pos;
static int8_t header_buf[SCF_HDR_LEN * SCF_TONES];
static int32_t last_preamble_weight;
static size_t data_bytes;
static size_t data_byte;
static size_t data_sym;
static int8_t data_buf[SCF_MSG_MAX][SCF_FEC_LEN * SCF_TONES];

static size_t cw_filter_idx;
static const float cw_filter_kernel[CW_FILTER_LEN] = {
    -0.058941394880239056,
    -0.067400917530924781,
    -0.074796392092324623,
    -0.080855873165259634,
    -0.085354127494019216,
    -0.088122662998354445,
    0.910942736322243651,
    -0.088122662998354445,
    -0.085354127494019216,
    -0.080855873165259634,
    -0.074796392092324623,
    -0.067400917530924781,
    -0.058941394880239056,
};

static bool initialized;

static float cw_filter(float *cw_filter_buf, float x)
{
    cw_filter_buf[cw_filter_idx] = x;

    float y = 0.0f;
    for (size_t i = 0; i < CW_FILTER_LEN; i++) {
        size_t pos = (cw_filter_idx + i + 1) % CW_FILTER_LEN;
        y += cw_filter_buf[pos] * cw_filter_kernel[i];
    }

    return y;
}

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

    float filtered_power[FFT_LEN];

    for (size_t i = 0; i < FFT_LEN; i++) {
        complex float bin = fft_buf[i];
        float power = crealf(bin) * crealf(bin) + cimagf(bin) * cimagf(bin);
        filtered_power[i] = cw_filter(&c->cw_filter_buf[i][0], power);
    }

    for (size_t i = 0; i < FFT_LEN; i++) {
        float filtered_power_sum = 0.0f;

        for (size_t j = 0; j < SCF_TONES; j++) {
            size_t f = (i + j * FFT_RATIO) % FFT_LEN;
            filtered_power_sum += fabs(filtered_power[f]);
        }
        for (size_t j = 0; j < SCF_TONES; j++) {
            size_t f = (i + j * FFT_RATIO) % FFT_LEN;
            c->demod_buf[i][demod_buf_idx][j] = 127.0f * filtered_power[f] / filtered_power_sum;
        }

        if (i == rx_bin) {
            int8_t max_weight = INT8_MIN;
            for (size_t j = 0; j < SCF_TONES; j++) {
                int8_t weight = c->demod_buf[i][demod_buf_idx][j];
                if (weight > max_weight) {
                    max_weight = weight;
                }
            }
        }
    }
}

static void find_preamble(struct sym_phase *p, int32_t *preamble_weight, size_t *preamble_bin)
{
    int32_t max_weight = INT32_MIN;
    size_t max_bin = 0;

    for (size_t b = 0; b < FFT_LEN; b++) {
        int32_t weight = 0;
        for (size_t i = 0; i < SCF_PREAMBLE; i++) {
            uint8_t t = preamble_ref[i];
            size_t pos = (demod_buf_idx + i + 1) % SCF_PREAMBLE;
            weight += p->demod_buf[b][pos][t];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b;
        }
    }

    *preamble_weight = max_weight;
    *preamble_bin = max_bin;
}

static bool decode_preamble(struct sym_phase *p, size_t preamble_bin)
{
    unsigned int correct = 0;

    for (size_t i = 0; i < SCF_PREAMBLE; i++) {
        size_t pos = (demod_buf_idx + i + 1) % SCF_PREAMBLE;
        int32_t max_weight = INT32_MIN;
        uint8_t symbol = 0;

        for (size_t t = 0; t < SCF_TONES; t++) {
            int32_t weight = p->demod_buf[preamble_bin][pos][t];
            if (weight > max_weight) {
                max_weight = weight;
                symbol = t;
            }
        }

        if (symbol == preamble_ref[i]) {
            correct++;
        }
    }

    if (correct >= SCF_PREAMBLE_REQUIRED) {
        return true;
    }

    return false;
}

static uint8_t ml_decode(int32_t *weight, int8_t *codeword, size_t codeword_size, uint8_t *code_table)
{
    int32_t max_weight = INT32_MIN;
    uint8_t max_symbol = 0;
    for (size_t i = 0; i < 256; i++) {
        int32_t weight = 0;
        for (size_t s = 0; s < codeword_size; s++) {
            size_t tone = code_table[i * codeword_size + s];
            weight += codeword[SCF_TONES * s + tone];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_symbol = i;
        }
    }

    *weight = max_weight;
    return max_symbol;
}

static void fft_window_init(void)
{
    for (size_t i = 0; i < BB_SYM_LEN; i++) {
        fft_window[i] = sin(M_PI * i / BB_SYM_LEN);
    }
}

void scf_rx_init(float freq, uint8_t *preamble)
{
    state = S_PREAMBLE;
    symbol_counter = 0;
    carrier_freq = freq;
    memcpy(preamble_ref, preamble, SCF_PREAMBLE);

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

    int32_t max_preamble_weight = INT32_MIN;
    size_t max_preamble_bin = 0;
    size_t max_preamble_phase = 0;

    for (size_t i = 0; i < SYM_PHASES; i++) {
        demodulate(&sym_phase[i]);

        int32_t preamble_weight;
        size_t preamble_bin;
        find_preamble(&sym_phase[i], &preamble_weight, &preamble_bin);

        if (preamble_weight > max_preamble_weight) {
            max_preamble_weight = preamble_weight;
            max_preamble_bin = preamble_bin;
            max_preamble_phase = i;
        }
    }

    bool preamble_found = decode_preamble(&sym_phase[max_preamble_phase], max_preamble_bin);

    switch (state) {
        case S_PREAMBLE:
            if (preamble_found) {
                last_preamble_weight = max_preamble_weight;
                header_pos = 0;
                rx_phase = max_preamble_phase;
                rx_bin = max_preamble_bin;
                state = S_HEADER;
            }

            break;

        case S_HEADER:
            if (header_pos == 0 && preamble_found && max_preamble_weight > last_preamble_weight) {
                last_preamble_weight = max_preamble_weight;
                rx_phase = max_preamble_phase;
                rx_bin = max_preamble_bin;
            } else {
                memcpy(&header_buf[header_pos * SCF_TONES], &sym_phase[rx_phase].demod_buf[rx_bin][demod_buf_idx][0], SCF_TONES);
                header_pos++;

                if (header_pos == SCF_HDR_LEN) {
                    int32_t header_weight;
                    uint8_t header_byte = ml_decode(&header_weight, header_buf, SCF_HDR_LEN, &scf_header_code[0][0]);
                    header_byte ^= 0xF0;

                    if (header_byte >= SCF_PACKET_TYPES) {
                        state = S_PREAMBLE;
                    } else {
                        data_bytes = scf_packet_len[header_byte];
                        data_byte = 0;
                        data_sym = 0;
                        state = S_DATA;
                    }
                }
            }

            break;

        case S_DATA:
            memcpy(&data_buf[data_byte][data_sym * SCF_TONES], &sym_phase[rx_phase].demod_buf[rx_bin][demod_buf_idx][0], SCF_TONES);

            if (data_sym == SCF_FEC_LEN - 1) {
                int32_t byte_weight;
                uint8_t byte_val = ml_decode(&byte_weight, &data_buf[data_byte][0], SCF_FEC_LEN, &scf_data_code[0][0]);
                msg[data_byte] = byte_val;
            }

            data_byte++;
            if (data_byte == data_bytes) {
                data_byte = 0;
                data_sym++;
            }

            if (data_sym == SCF_FEC_LEN) {
                msg_len = data_bytes;
                state = S_PREAMBLE;
            }

            break;
    }


    demod_buf_idx = (demod_buf_idx + 1) % SCF_PREAMBLE;
    cw_filter_idx = (cw_filter_idx + 1) % CW_FILTER_LEN;

    symbol_counter++;

    return msg_len;
}
