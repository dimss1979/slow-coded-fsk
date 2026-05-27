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
    complex float *source;
    uint32_t decoded_symbol[SCF_MAX_FEC_LEN];
    size_t prev_corr_peak_index;
};

float center_freq;
float phase;
complex float baseband[SCF_SYM_LEN * 2];
complex float fir_tail[SCF_FIR_LEN_RF];

static complex float zadoff_chu_template[SCF_SYM_LEN];
static struct sym_phase sym_phase[SYM_PHASES];
static fftwf_plan fft_plan, ifft_plan;
static fftwf_complex *fft_buf;
static unsigned int symbol_counter;
static unsigned int symbol_skip;
static size_t decoded_symbol_index;

static bool initialized;

static void downconvert(float *signal)
{
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        baseband[i] = baseband[i + SCF_SYM_LEN];
    }

    complex float baseband_mixed[SCF_SYM_LEN] = {0};

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        float carrier_i = sinf(phase);
        float carrier_q = cosf(phase);
        baseband_mixed[i] = signal[i] * carrier_i + signal[i] * I * carrier_q;

        phase += 2.0f * M_PI * center_freq * (1.0f / (float) SCF_SRATE);
        while (phase > 2.0f * M_PI) {
            phase -= 2.0f * M_PI;
        }
    }
    scf_filter_rf(baseband + SCF_SYM_LEN, baseband_mixed, fir_tail);
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
    complex float correlated_td[SCF_SYM_LEN] = {0};

    correlate_against_template(correlated_td, c->source, zadoff_chu_template);

    float max_peak_value = 0.0f;
    size_t max_peak_index = 0;
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        complex float peak = correlated_td[i];
        float peak_value = peak * conjf(peak);
        if (peak_value > max_peak_value) {
            max_peak_value = peak_value;
            max_peak_index = i;
        }
    }

    float max_symbol_peak_value = 0.0f;
    size_t max_symbol_index = 0;
    for (uint32_t symbol = 0; symbol < 256; symbol++) {
        size_t symbol_offset = symbol * SCF_DEC_RATIO;
        size_t peak_index = (c->prev_corr_peak_index - symbol_offset + SCF_SYM_LEN) % SCF_SYM_LEN;
        complex float peak = correlated_td[peak_index];
        float peak_value = peak * conjf(peak);
        if (peak_value > max_symbol_peak_value) {
            max_symbol_peak_value = peak_value;
            max_symbol_index = symbol;
        }
    }

    c->decoded_symbol[decoded_symbol_index] = max_symbol_index;
    c->prev_corr_peak_index = max_peak_index;
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
        sym_phase[i].source = &baseband[i * SCF_SYM_LEN / SYM_PHASES];
        // To prevent Reed Solomon decoder from decoding the symbol prematurely
        getrandom(sym_phase[i].decoded_symbol, sizeof(sym_phase[i].decoded_symbol), 0);
    }

    center_freq = freq;
    phase = 0.0f;

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

    downconvert(signal);

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
