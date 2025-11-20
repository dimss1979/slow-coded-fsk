#include <complex.h>
#include <math.h>
#include <stddef.h>

#include "scf.h"
#include "scf_inner.h"
#include "scf_filter.h"
#include "scf_tx.h"

#define MOD_FILTER_LEN 8

static complex float fir_tail[SCF_FIR_LEN_RF];
static float carrier_freq;
static float carrier_phase;
static float bb_phase;
static const float freq_step = (float) SCF_SRATE / SCF_CHIP_LEN;

static float mod_filter_buf[MOD_FILTER_LEN];
static const float mod_filter_kernel[MOD_FILTER_LEN] = {
    0.012564432734645786,
    0.058060308288456049,
    0.161078275758929301,
    0.268296983217968910,
    0.268296983217968910,
    0.161078275758929301,
    0.058060308288456049,
    0.012564432734645786,
};
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

void scf_tx_init(float freq)
{
    scf_filter_init();

    carrier_freq = freq;
}

void scf_tx(float *passband, uint8_t symbol, float gain)
{
    size_t passband_idx = 0;
    for (size_t c = 0; c < SCF_CHIPS; c++) {
        float freq = freq_step * scf_inner_code[symbol][c];
        complex float baseband[SCF_CHIP_LEN] = {0};
        complex float filtered[SCF_CHIP_LEN];

        for (size_t i = 0; i < SCF_BB_CHIP_LEN; i++) {
            baseband[i * SCF_DEC_RATIO] = sinf(bb_phase) + I * cosf(bb_phase);

            bb_phase += 2.0f * M_PI * mod_filter(freq) * (1.0f / SCF_BB_SRATE);
            while (bb_phase > 2.0f * M_PI) {
                bb_phase -= 2.0f * M_PI;
            }
            while (bb_phase < 2.0f * M_PI) {
                bb_phase += 2.0f * M_PI;
            }
        }

        scf_filter_rf(filtered, baseband, fir_tail);

        for (size_t i = 0; i < SCF_CHIP_LEN; i++) {
            passband[passband_idx] = sinf(carrier_phase) * crealf(filtered[i]) - cosf(carrier_phase) * cimagf(filtered[i]);
            passband[passband_idx] *= gain;
            passband_idx++;

            carrier_phase += 2.0f * M_PI * carrier_freq * (1.0f / SCF_SRATE);
            while (carrier_phase > 2.0f * M_PI) {
                carrier_phase -= 2.0f * M_PI;
            }
        }
    }
}
