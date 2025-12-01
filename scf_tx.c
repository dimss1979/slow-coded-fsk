#include <complex.h>
#include <math.h>
#include <stddef.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_tx.h"

#define MOD_FILTER_LEN (SCF_SYM_LEN / 10)

static float carrier_freq;
static float carrier_phase;
static const float freq_step = (float) SCF_SRATE / SCF_SYM_LEN;

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

void mod_filter_init(void)
{
    float dc_gain = 0.0f;

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        mod_filter_kernel[i] = sin(M_PI * i / MOD_FILTER_LEN);
        dc_gain += mod_filter_kernel[i];
    }

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        mod_filter_kernel[i] /= dc_gain;
    }
}

void scf_tx_init(float freq)
{
    scf_filter_init();
    mod_filter_init();

    carrier_freq = freq;
}

void scf_tx(float *passband, uint8_t symbol, float gain)
{
    float freq = carrier_freq + freq_step * symbol;

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        passband[i] = gain * sinf(carrier_phase);

        carrier_phase += 2.0f * M_PI * mod_filter(freq) * (1.0f / SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }
    }
}
