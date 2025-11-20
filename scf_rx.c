#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "scf.h"
#include "scf_inner.h"
#include "scf_filter.h"
#include "scf_rx.h"

#define CHIP_PHASES 8
#define CHIP_FREQS 4
#define FREQ_STEP_PER_BIN 4
#define CW_FILTER_LEN 13

#define TRACK_LEN (9)

struct rx_chain {
    float cw_filter_buf[SCF_BB_CHIP_LEN][CW_FILTER_LEN];
    int8_t mfsk_history[SCF_BB_CHIP_LEN][SCF_CHIPS * SCF_PREAMBLE][SCF_FREQS];
    int32_t preamble_weight;
    size_t preamble_bin;
    int32_t weight_history[TRACK_LEN][SCF_CHIPS];
    struct scf_soft_symbol symbol_history[SCF_CHIPS];
    struct scf_soft_symbol decoded_symbol;
    int32_t decoded_phase_weight;
    size_t decoded_phase;
};

struct rx_worker {
    complex float *source;
    fftwf_complex *sa_fft_buf;
    struct rx_chain rx_chains[CHIP_FREQS];
    struct scf_soft_symbol *decoded_symbol;
};

static struct rx_worker *rx_workers;
static fftwf_plan sa_fft;
static complex float input_chips[SCF_BB_CHIP_LEN * 2];
static float carrier_freq;
static float carrier_phase;
static complex float fir_tail[SCF_FIR_LEN_RF];
static complex float shifter_wavetable[CHIP_FREQS][SCF_BB_CHIP_LEN];
static uint8_t rx_preamble[SCF_PREAMBLE];
static size_t rx_bin;
static size_t mfsk_idx;
static size_t inner_idx;
static size_t outer_idx;
static size_t chip_cnt;

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

static void shift_input_chips(float *chip)
{
    memcpy(&input_chips[0], &input_chips[SCF_BB_CHIP_LEN], SCF_BB_CHIP_LEN * sizeof(input_chips[0]));

    complex float baseband[SCF_CHIP_LEN];
    complex float baseband_filtered[SCF_CHIP_LEN];
    for (size_t i = 0; i < SCF_CHIP_LEN; i++) {
        float carrier_i = sinf(carrier_phase);
        float carrier_q = cosf(carrier_phase);
        baseband[i] = chip[i] * carrier_i + chip[i] * I * carrier_q;

        carrier_phase += 2.0f * M_PI * carrier_freq * (1.0f / (float) SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }
    }
    scf_filter_rf(baseband_filtered, baseband, fir_tail);
    for (size_t i = 0; i < SCF_BB_CHIP_LEN; i++) {
        input_chips[SCF_BB_CHIP_LEN + i] = baseband_filtered[i * SCF_DEC_RATIO];
    }
}

static void update_mfsk_history(struct rx_worker *w, struct rx_chain *c, size_t freq_i)
{
    for (size_t i = 0; i < SCF_BB_CHIP_LEN; i++) {
        w->sa_fft_buf[i] = w->source[i] * shifter_wavetable[freq_i][i];
    }

    fftwf_execute_dft(sa_fft, w->sa_fft_buf, w->sa_fft_buf);

    float filtered_power[SCF_BB_CHIP_LEN];

    for (size_t i = 0; i < SCF_BB_CHIP_LEN; i++) {
        complex float bin = w->sa_fft_buf[i];
        float power = crealf(bin) * crealf(bin) + cimagf(bin) * cimagf(bin);
        filtered_power[i] = cw_filter(&c->cw_filter_buf[i][0], power);
    }

    for (size_t i = 0; i < SCF_BB_CHIP_LEN; i++) {
        float filtered_power_sum = 0.0f;

        for (size_t j = 0; j < SCF_FREQS; j++) {
            size_t f = (i + j) % SCF_BB_CHIP_LEN;
            filtered_power_sum += fabs(filtered_power[f]);
        }
        for (size_t j = 0; j < SCF_FREQS; j++) {
            size_t f = (i + j) % SCF_BB_CHIP_LEN;
            c->mfsk_history[i][mfsk_idx][j] = 127.0f * filtered_power[f] / filtered_power_sum;
        }
    }
}

static void update_preamble(struct rx_chain *c)
{
    int32_t max_weight = INT32_MIN;
    size_t max_bin = 0;

    for (size_t b = 0; b < SCF_BB_CHIP_LEN; b++) {
        int32_t weight = 0;
        size_t preamble_pos = mfsk_idx;
        for (size_t i = 0; i < SCF_PREAMBLE; i++) {
            size_t s = rx_preamble[i];
            for (size_t t = 0; t < SCF_CHIPS; t++) {
                preamble_pos = (preamble_pos + 1) % (SCF_CHIPS * SCF_PREAMBLE);
                weight += c->mfsk_history[b][preamble_pos][scf_inner_code[s][t]];
            }
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_bin = b;
        }
    }

    c->preamble_weight = max_weight;
    c->preamble_bin = max_bin;
}

static void decode_preamble(struct rx_chain *c)
{
    unsigned int correct = 0;
    int32_t preamble_weight[SCF_PREAMBLE];

    for (size_t i = 0; i < SCF_PREAMBLE; i++) {
        int32_t max_weight = INT32_MIN;
        uint8_t max_symbol = 0;

        for (size_t s = 0; s < SCF_SYMBOL_M; s++) {
            int32_t weight = 0;

            for (size_t t = 0; t < SCF_CHIPS; t++) {
                size_t preamble_pos = (mfsk_idx + 1 + (i * SCF_CHIPS) + t) % (SCF_CHIPS * SCF_PREAMBLE);
                weight += c->mfsk_history[c->preamble_bin][preamble_pos][scf_inner_code[s][t]];
            }

            if (weight > max_weight) {
                max_weight = weight;
                max_symbol = s;
            }
        }

        preamble_weight[i] = max_weight;
        if (max_symbol == rx_preamble[i]) {
            correct++;
        }
    }

    if (correct >= SCF_PREAMBLE_REQUIRED) {
        rx_bin = c->preamble_bin;
        chip_cnt = SCF_CHIPS * 3 / 2;

        for (size_t i = 0; i < SCF_PREAMBLE; i++) {
            size_t history_pos = (outer_idx - SCF_PREAMBLE + 1 + i + TRACK_LEN) % TRACK_LEN;
            c->weight_history[history_pos][inner_idx] = preamble_weight[i];
        }
    }
}

static void decode_current_symbol(struct rx_chain *c)
{
    int32_t max_weight = INT32_MIN;
    uint8_t max_symbol = 0;

    for (size_t s = 0; s < SCF_SYMBOL_M; s++) {
        int32_t weight = 0;
        for (size_t t = 0; t < SCF_CHIPS; t++) {
            size_t history_pos = (mfsk_idx + 1 + t + (SCF_CHIPS * (SCF_PREAMBLE - 1))) % (SCF_CHIPS * SCF_PREAMBLE);
            weight += c->mfsk_history[rx_bin][history_pos][scf_inner_code[s][t]];
        }

        if (weight > max_weight) {
            max_weight = weight;
            max_symbol = s;
        }
    }

    c->weight_history[outer_idx][inner_idx] = max_weight;
    c->symbol_history[inner_idx].symbol = max_symbol;
    c->symbol_history[inner_idx].weight = max_weight;

}

static void find_peak_symbol(struct rx_chain *c)
{
    int32_t max_phase_weight = INT32_MIN;
    size_t max_i = 0;

    for (size_t i = 0; i < SCF_CHIPS; i++) {
        int32_t phase_weight = 0;
        for (size_t j = 0; j < TRACK_LEN; j++) {
            phase_weight += c->weight_history[j][i];
        }

        if (phase_weight > max_phase_weight) {
            max_phase_weight = phase_weight;
            max_i = i;
        }
    }
    c->decoded_symbol = c->symbol_history[max_i];
    c->decoded_phase_weight = max_phase_weight;

    size_t phase = (inner_idx - max_i + SCF_CHIPS) % SCF_CHIPS;
    c->decoded_phase = phase;
}

static void rx_worker_preamble(struct rx_worker *w)
{
    for (size_t i = 0; i < CHIP_FREQS; i++) {
        struct rx_chain *c = &w->rx_chains[i];
        update_mfsk_history(w, c, i);
        update_preamble(c);
    }
}

static void rx_worker_data(struct rx_worker *w)
{
    for (size_t i = 0; i < CHIP_FREQS; i++) {
        struct rx_chain *c = &w->rx_chains[i];
        decode_current_symbol(c);
        if (chip_cnt == 1) {
            find_peak_symbol(c);
        }
    }
}

static void shifter_wavetable_init(void)
{
    float freq_step = 1.0f / ((float) FREQ_STEP_PER_BIN * SCF_CHIP_LEN / SCF_SRATE);

    for (int freq_i = 0; freq_i < CHIP_FREQS; freq_i++) {
        float shift_freq = freq_step * (freq_i - CHIP_FREQS / 2);
        float shift_phase = 0.0f;
        for (size_t i = 0; i < SCF_BB_CHIP_LEN; i++) {
            shifter_wavetable[freq_i][i] = sinf(shift_phase) + cosf(shift_phase) * I;
            shift_phase += 2.0f * M_PI * shift_freq * (1.0f / SCF_BB_SRATE);
            while (shift_phase > 2.0f * M_PI) {
                shift_phase -= 2.0f * M_PI;
            }
            while (shift_phase < -2.0f * M_PI) {
                shift_phase += 2.0f * M_PI;
            }
        }
    }
}

void scf_rx_init(float freq, uint8_t *preamble)
{
    carrier_freq = freq;
    memcpy(rx_preamble, preamble, SCF_PREAMBLE);

    if (!rx_workers) {
        rx_workers = calloc(CHIP_PHASES, sizeof(*rx_workers));
        assert(rx_workers);

        for (size_t i = 0; i < CHIP_PHASES; i++) {
            struct rx_worker *w = &rx_workers[i];

            w->source = &input_chips[SCF_BB_CHIP_LEN - i * SCF_BB_CHIP_LEN / CHIP_PHASES];
            w->sa_fft_buf = fftwf_alloc_complex(SCF_BB_CHIP_LEN);
            assert(w->sa_fft_buf);
        }

        sa_fft = fftwf_plan_dft_1d(
            SCF_BB_CHIP_LEN,
            rx_workers[0].sa_fft_buf,
            rx_workers[0].sa_fft_buf,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        assert(sa_fft);

        mfsk_idx = 0;
        inner_idx = 0;
        outer_idx = 0;
        chip_cnt = SCF_CHIPS;
        shifter_wavetable_init();
    }
}

bool scf_rx_chip(struct scf_soft_symbol *symbol, float *chip)
{
    shift_input_chips(chip);

    for (size_t i = 0; i < CHIP_PHASES; i++) {
        struct rx_worker *w = &rx_workers[i];
        rx_worker_preamble(w);
    }

    int32_t max_preamble_weight = INT32_MIN;
    struct rx_chain *preamble_chain = NULL;

    for (size_t i = 0; i < CHIP_PHASES; i++) {
        struct rx_worker *w = &rx_workers[i];
        for (size_t i = 0; i < CHIP_FREQS; i++) {
            struct rx_chain *c = &w->rx_chains[i];
            if (c->preamble_weight > max_preamble_weight) {
                max_preamble_weight = c->preamble_weight;
                preamble_chain = c;
            }
        }
    }

    decode_preamble(preamble_chain);

    for (size_t i = 0; i < CHIP_PHASES; i++) {
        struct rx_worker *w = &rx_workers[i];
        rx_worker_data(w);
    }

    mfsk_idx = (mfsk_idx + 1) % (SCF_CHIPS * SCF_PREAMBLE);
    inner_idx = (inner_idx + 1) % SCF_CHIPS;
    if (inner_idx == 0) {
        outer_idx = (outer_idx + 1) % TRACK_LEN;
    }
    cw_filter_idx = (cw_filter_idx + 1) % CW_FILTER_LEN;
    chip_cnt--;

    if (!chip_cnt) {
        int32_t max_phase_weight = 0;
        size_t best_phase = 0;

        for (size_t i = 0; i < CHIP_PHASES; i++) {
            struct rx_worker *w = &rx_workers[i];
            for (size_t i = 0; i < CHIP_FREQS; i++) {
                struct rx_chain *c = &w->rx_chains[i];
                if (c->decoded_phase_weight > max_phase_weight) {
                    max_phase_weight = c->decoded_phase_weight;
                    *symbol = c->decoded_symbol;
                    best_phase = c->decoded_phase;
                }
            }
        }

        if (best_phase < SCF_CHIPS / 2) {
            chip_cnt = SCF_CHIPS + 1;
        } else if (best_phase > SCF_CHIPS / 2) {
            chip_cnt = SCF_CHIPS - 1;
        } else {
            chip_cnt = SCF_CHIPS;
        }
        return true;
    } else {
        return false;
    }
}
