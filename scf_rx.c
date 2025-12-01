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
    int8_t preamble_buf[FFT_LEN][SCF_PREAMBLE][SCF_TONES];
    int32_t preamble_weight;
    size_t preamble_bin;
};

struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan;
fftwf_complex *fft_buf;
static complex float input_signal[BB_SYM_LEN * 2];
static float carrier_freq;
static float carrier_phase;
static complex float fir_tail[SCF_FIR_LEN_RF];
static float fft_window[BB_SYM_LEN];
static uint8_t preamble_ref[SCF_PREAMBLE];
static size_t preamble_idx;
static size_t signal_bin;

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

bool initialized;

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
            c->preamble_buf[i][preamble_idx][j] = 127.0f * filtered_power[f] / filtered_power_sum;
        }
    }
}

static void find_preamble(struct sym_phase *p)
{
    int32_t max_weight = INT32_MIN;
    size_t max_bin = 0;

    for (size_t b = 0; b < FFT_LEN; b++) {
        int32_t weight = 0;
        for (size_t i = 0; i < SCF_PREAMBLE; i++) {
            uint8_t t = preamble_ref[i];
            size_t pos = (preamble_idx + i + 1) % SCF_PREAMBLE;
            weight += p->preamble_buf[b][pos][t];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b;
        }
    }

    p->preamble_weight = max_weight;
    p->preamble_bin = max_bin;
}

static void decode_preamble(struct sym_phase *p)
{
    unsigned int correct = 0;

    for (size_t i = 0; i < SCF_PREAMBLE; i++) {
        size_t pos = (preamble_idx + i + 1) % SCF_PREAMBLE;
        int32_t max_weight = INT32_MIN;
        uint8_t symbol = 0;

        for (size_t t = 0; t < SCF_TONES; t++) {
            int32_t weight = p->preamble_buf[p->preamble_bin][pos][t];
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
        printf("Preamble decoded\n");
        signal_bin = p->preamble_bin;
    }
}

static void fft_window_init(void)
{
    for (size_t i = 0; i < BB_SYM_LEN; i++) {
        fft_window[i] = sin(M_PI * i / BB_SYM_LEN);
    }
}

void scf_rx_init(float freq, uint8_t *preamble)
{
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

    for (size_t i = 0; i < SYM_PHASES; i++) {
        demodulate(&sym_phase[i]);
        find_preamble(&sym_phase[i]);
    }

    int32_t max_preamble_weight = INT32_MIN;
    struct sym_phase *preamble_phase = NULL;

    for (size_t i = 0; i < SYM_PHASES; i++) {
        struct sym_phase *p = &sym_phase[i];
        if (p->preamble_weight > max_preamble_weight) {
            max_preamble_weight = p->preamble_weight;
            preamble_phase = p;
        }
    }

    decode_preamble(preamble_phase);

    preamble_idx = (preamble_idx + 1) % SCF_PREAMBLE;
    cw_filter_idx = (cw_filter_idx + 1) % CW_FILTER_LEN;

    return msg_len;
}
