#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include "scf_private.h"

#define MOD_FILTER_LEN (SCF_BB_SYM_LEN / 20)
#define FREQ_STEP ((float)SCF_BB_SRATE / SCF_BB_SYM_LEN)

static float carrier_freq;
static float carrier_phase;
static float baseband_phase;

static float mod_filter_buf[MOD_FILTER_LEN];
static float mod_filter_kernel[MOD_FILTER_LEN];
static size_t mod_filter_idx;

static float mod_filter(float x)
{
    mod_filter_buf[mod_filter_idx] = x;
    mod_filter_idx = (mod_filter_idx + 1) % MOD_FILTER_LEN;

    float y = 0.0f;
    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        size_t pos = (mod_filter_idx + i) % MOD_FILTER_LEN;
        y += mod_filter_buf[pos] * mod_filter_kernel[i];
    }

    return y;
}

static void mod_filter_init(void)
{
    float dc_gain = 0.0f;

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        int mod_filter_len_i = MOD_FILTER_LEN;
        float mod_filter_len_f = (float) mod_filter_len_i;
        mod_filter_kernel[i] = sinf((M_PI * i) / mod_filter_len_f);
        dc_gain += mod_filter_kernel[i];
    }

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        mod_filter_kernel[i] /= dc_gain;
    }
}

void scf_tx_init(float freq)
{
    mod_filter_init();
    memset(mod_filter_buf, 0, sizeof(mod_filter_buf));
    mod_filter_idx = 0;
    scf_filter_reset_tx();

    carrier_freq = freq;
    carrier_phase = 0.0f;
    baseband_phase = 0.0f;
}

void scf_tx(float signal[SCF_SYM_LEN], uint32_t symbol, float gain)
{
    complex float baseband[SCF_SYM_LEN] = {0};
    complex float baseband_filtered[SCF_SYM_LEN] = {0};

    float freq = FREQ_STEP * ((float) symbol + 0.5f - (float) SCF_TONES / 2.0f);

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        baseband[i * SCF_DEC_RATIO] = gain * (complex float) (sinf(baseband_phase) + I * cosf(baseband_phase));

        baseband_phase += 2.0f * M_PI * mod_filter(freq) * (1.0f / SCF_BB_SRATE);
        while (baseband_phase > 2.0f * M_PI) {
            baseband_phase -= 2.0f * M_PI;
        }
        while (baseband_phase < -2.0f * M_PI) {
            baseband_phase += 2.0f * M_PI;
        }
    }

    scf_filter_tx(baseband_filtered, baseband);

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        carrier_phase += 2.0f * M_PI * carrier_freq * (1.0f / SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }

        complex float bb = baseband_filtered[i];
        signal[i] = crealf(bb) * sinf(carrier_phase) - cimagf(bb) * cosf(carrier_phase);
    }
}
