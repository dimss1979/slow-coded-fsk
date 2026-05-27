#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <sys/random.h>

#include "scf.h"
#include "scf_private.h"

static float center_freq;
static float phase;
static size_t code_shift;

static complex float fir_tail[SCF_FIR_LEN_RF];

void scf_tx_init(float freq)
{
    center_freq = freq;
}

void scf_tx(float signal[SCF_SYM_LEN], uint32_t symbol, float gain)
{
    code_shift = (code_shift + symbol) % SCF_BB_SYM_LEN;

    complex float baseband[SCF_SYM_LEN] = {0};
    complex float baseband_filtered[SCF_SYM_LEN] = {0};

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        size_t code_i = (i + code_shift) % SCF_BB_SYM_LEN;
        baseband[i * SCF_DEC_RATIO] = scf_packet_zadoff_chu_sequence[code_i];
    }

    scf_filter_rf(baseband_filtered, baseband, fir_tail);

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        phase += 2.0f * M_PI * center_freq * (1.0f / SCF_SRATE);
        while (phase > 2.0f * M_PI) {
            phase -= 2.0f * M_PI;
        }

        complex float bb = baseband_filtered[i];
        signal[i] = crealf(bb) * sinf(phase) + cimagf(bb) * cosf(phase);

        signal[i] *= gain;
    }
}
