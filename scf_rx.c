#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include "scf_private.h"

#define SYM_PHASES 8

struct sym_phase {
    complex float *source;
};

static struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan;
static fftwf_complex *fft_buf;
static complex float input_signal[SCF_BB_SYM_LEN * 2];
static float carrier_freq;
static float carrier_phase;
static complex float fir_tail[SCF_FIR_LEN_RF];
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

static void demodulate(struct sym_phase *c)
{
    // TODO: demodulate symbol
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

        fft_buf = fftwf_alloc_complex(SCF_SYM_LEN);
        assert(fft_buf);

        fft_plan = fftwf_plan_dft_1d(
            SCF_SYM_LEN,
            fft_buf,
            fft_buf,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        assert(fft_plan);

        initialized = true;
    }
}

void scf_rx(scf_rx_result *result, float signal[SCF_SYM_LEN])
{
    memset(result, 0, sizeof(scf_rx_result));
    result->symbol_counter = symbol_counter;

    downconvert(signal);

    for (size_t i = 0; i < SYM_PHASES; i++) {
        demodulate(&sym_phase[i]);

        for (size_t pt = 0; pt < SCF_PACKET_TYPES; pt++) {
            // TODO: decode packet

            uint8_t outer_codeword[SCF_MAX_FEC_LEN] = {0};
            memset(outer_codeword, 123, sizeof(outer_codeword));
            if (!result->got_msg && scf_packet_decode(result, outer_codeword, pt)) {
                symbol_skip = 5;
                break;
            }
        }
    }

    symbol_counter++;
    if (symbol_skip) {
        symbol_skip--;
    }
}
