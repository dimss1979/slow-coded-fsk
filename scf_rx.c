#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <sys/random.h>

#include "scf.h"
#include "scf_private.h"

#define SYM_PHASES 1

struct sym_phase {
    complex float *pilot_source;
    complex float *bearer_source;
    uint32_t decoded_symbol[SCF_MAX_FEC_LEN];
};

struct signal_chain {
    float freq;
    float phase;
    complex float baseband[SCF_SYM_LEN * 2];
    complex float fir_tail[SCF_FIR_LEN_RF];
} pilot_chain, bearer_chain;

static complex float zadoff_chu_template[SCF_SYM_LEN];
static struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan, ifft_plan;
static fftwf_complex *fft_buf;
static unsigned int symbol_counter;
static unsigned int symbol_skip;
static size_t decoded_symbol_index;

static bool initialized;

static void downconvert(struct signal_chain *chain, float *signal)
{
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        chain->baseband[i] = chain->baseband[i + SCF_SYM_LEN];
    }

    complex float baseband[SCF_SYM_LEN];
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        float carrier_i = sinf(chain->phase);
        float carrier_q = cosf(chain->phase);
        baseband[i] = signal[i] * carrier_i + signal[i] * I * carrier_q;

        chain->phase += 2.0f * M_PI * chain->freq * (1.0f / (float) SCF_SRATE);
        while (chain->phase > 2.0f * M_PI) {
            chain->phase -= 2.0f * M_PI;
        }
    }
    scf_filter_rf(chain->baseband + SCF_SYM_LEN, baseband, chain->fir_tail);
}

static void correlate_against_template(complex float *out_td, complex float *in_td, complex float *template_fd)
{
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        fft_buf[i] = in_td[i];
    }
    fftwf_execute(fft_plan);
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        fft_buf[i] = fft_buf[i] * template_fd[i];
    }
    fftwf_execute(ifft_plan);
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        out_td[i] = fft_buf[i] / (float) SCF_SYM_LEN;
    }
}

static void demodulate(struct sym_phase *c)
{
    complex float correlated_pilot_td[SCF_SYM_LEN] = {0};
    complex float correlated_bearer_td[SCF_SYM_LEN] = {0};

    correlate_against_template(correlated_pilot_td, c->pilot_source, zadoff_chu_template);
    correlate_against_template(correlated_bearer_td, c->bearer_source, zadoff_chu_template);

    float max_pilot_peak_value = 0.0f;
    size_t max_pilot_peak_index = 0;
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        complex float peak = correlated_pilot_td[i];
        float peak_value = peak * conjf(peak);
        if (peak_value > max_pilot_peak_value) {
            max_pilot_peak_value = peak_value;
            max_pilot_peak_index = i;
        }
    }


    float max_bearer_peak_value = 0.0f;
    size_t max_bearer_peak_symbol_index = 0;
    for (uint32_t symbol = 0; symbol < 256; symbol++) {
        size_t bearer_offset = (symbol + 1) * SCF_DEC_RATIO;
        size_t bearer_index = (max_pilot_peak_index - bearer_offset+ SCF_SYM_LEN) % SCF_SYM_LEN;
        complex float bearer_peak = correlated_bearer_td[bearer_index];
        float bearer_peak_value = bearer_peak * conjf(bearer_peak);
        if (bearer_peak_value > max_bearer_peak_value) {
            max_bearer_peak_value = bearer_peak_value;
            max_bearer_peak_symbol_index = symbol;
        }
    }

    c->decoded_symbol[decoded_symbol_index] = max_bearer_peak_symbol_index;
}

static void generate_zadoff_chu_template(void)
{
    complex float fir_tail[SCF_FIR_LEN_RF] = {0};
    complex float upconverted[SCF_SYM_LEN] = {0};

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        upconverted[i * SCF_DEC_RATIO] = scf_packet_zadoff_chu_sequence[i];
    }

    scf_filter_rf(fft_buf, upconverted, fir_tail);
    for (size_t i = 0; i < SCF_FIR_LEN_RF; i++) {
        fft_buf[i] += fir_tail[i];
    }

    fftwf_execute(fft_plan);
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        zadoff_chu_template[i] = conjf(fft_buf[i]);
    }
}

void scf_rx_init(float freq)
{
    symbol_counter = 0;
    symbol_skip = 0;

    for (size_t i = 0; i < SYM_PHASES; i++) {
        sym_phase[i].pilot_source = &pilot_chain.baseband[i * SCF_SYM_LEN / SYM_PHASES];
        sym_phase[i].bearer_source = &bearer_chain.baseband[i * SCF_SYM_LEN / SYM_PHASES];
        // To prevent Reed Solomon decoder from decoding the symbol prematurely
        getrandom(sym_phase[i].decoded_symbol, sizeof(sym_phase[i].decoded_symbol), 0);
    }

    memset(&pilot_chain, 0, sizeof(pilot_chain));
    memset(&bearer_chain, 0, sizeof(bearer_chain));
    pilot_chain.freq = freq - SCF_GUARD_BAND / 2 - SCF_BB_SRATE / 2;
    bearer_chain.freq = freq + SCF_GUARD_BAND / 2 + SCF_BB_SRATE / 2;

    if (!initialized) {
        fft_buf = fftwf_alloc_complex(SCF_SYM_LEN);
        assert(fft_buf);

        fft_plan = fftwf_plan_dft_1d(
            SCF_SYM_LEN,
            fft_buf,
            fft_buf,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        ifft_plan = fftwf_plan_dft_1d(
            SCF_SYM_LEN,
            fft_buf,
            fft_buf,
            FFTW_BACKWARD,
            FFTW_ESTIMATE
        );
        assert(fft_plan);
        assert(ifft_plan);

        initialized = true;
    }

    generate_zadoff_chu_template();
}

void scf_rx(scf_rx_result *result, float signal[SCF_SYM_LEN])
{
    memset(result, 0, sizeof(scf_rx_result));
    result->symbol_counter = symbol_counter;

    downconvert(&pilot_chain, signal);
    downconvert(&bearer_chain, signal);

    for (size_t phase_index = 0; phase_index < SYM_PHASES; phase_index++) {
        struct sym_phase *c = &sym_phase[phase_index];
        demodulate(c);

        if (symbol_skip) {
            continue;
        }

        for (size_t pt = 0; pt < SCF_PACKET_TYPES; pt++) {
            size_t packet_len = scf_packet_fec_len[pt];
            uint8_t packet[packet_len];
            for (size_t i = 0; i < packet_len; i++) {
                size_t pos = (i + decoded_symbol_index - packet_len + 1 + SCF_MAX_FEC_LEN) % SCF_MAX_FEC_LEN;
                packet[i] = c->decoded_symbol[pos];
            }

            if (!result->got_msg && scf_packet_decode(result, packet, pt)) {
                symbol_skip = 5;
                break;
            }
        }
    }

    symbol_counter++;
    decoded_symbol_index = (decoded_symbol_index + 1) % SCF_MAX_FEC_LEN;

    if (symbol_skip) {
        symbol_skip--;
    }
}
