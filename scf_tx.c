#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <sys/random.h>
#include "scf.h"
#include "scf_private.h"

static float pilot_freq;
static float bearer_freq;
static float pilot_phase;
static float bearer_phase;

static complex float pilot_fir_tail[SCF_FIR_LEN_RF];
static complex float bearer_fir_tail[SCF_FIR_LEN_RF];

void scf_tx_init(float freq)
{
    pilot_freq = freq - SCF_GUARD_BAND / 2 - SCF_BB_SRATE / 2;
    bearer_freq = freq + SCF_GUARD_BAND / 2 + SCF_BB_SRATE / 2;
}

void scf_tx(float signal[SCF_SYM_LEN], uint32_t symbol, float gain)
{
    complex float pilot_baseband[SCF_SYM_LEN] = {0};
    complex float pilot_baseband_filtered[SCF_SYM_LEN] = {0};
    complex float bearer_baseband[SCF_SYM_LEN] = {0};
    complex float bearer_baseband_filtered[SCF_SYM_LEN] = {0};

    uint32_t random_number = 0;
    getrandom(&random_number, sizeof(random_number), 0);
    size_t pilot_offset = random_number % SCF_BB_SYM_LEN;
    size_t bearer_offset = (pilot_offset + symbol + 1);

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        size_t pilot_i = (i + pilot_offset) % SCF_BB_SYM_LEN;
        size_t bearer_i = (i + bearer_offset) % SCF_BB_SYM_LEN;
        pilot_baseband[i * SCF_DEC_RATIO] = scf_packet_zadoff_chu_sequence[pilot_i];
        bearer_baseband[i * SCF_DEC_RATIO] = scf_packet_zadoff_chu_sequence[bearer_i];
    }

    scf_filter_rf(pilot_baseband_filtered, pilot_baseband, pilot_fir_tail);
    scf_filter_rf(bearer_baseband_filtered, bearer_baseband, bearer_fir_tail);

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        pilot_phase += 2.0f * M_PI * pilot_freq * (1.0f / SCF_SRATE);
        while (pilot_phase > 2.0f * M_PI) {
            pilot_phase -= 2.0f * M_PI;
        }

        bearer_phase += 2.0f * M_PI * bearer_freq * (1.0f / SCF_SRATE);
        while (bearer_phase > 2.0f * M_PI) {
            bearer_phase -= 2.0f * M_PI;
        }

        complex float pilot_bb = pilot_baseband_filtered[i];
        complex float bearer_bb = bearer_baseband_filtered[i];
        signal[i] = crealf(pilot_bb) * sinf(pilot_phase) - cimagf(pilot_bb) * cosf(pilot_phase) +
                    crealf(bearer_bb) * sinf(bearer_phase) - cimagf(bearer_bb) * cosf(bearer_phase);

        signal[i] *= gain;
    }
}
