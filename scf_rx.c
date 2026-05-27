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

#define SYM_PHASES 8
#define CFO_STEPS  4

// min (257+256)*8 = 4104
#define FFT_LEN    4320

#define RX_CHAINS (SYM_PHASES * CFO_STEPS)

struct rx_chain {
    size_t cfo_index;
    complex float *source;
    uint32_t decoded_symbol[SCF_MAX_FEC_LEN];
    size_t prev_corr_peak_index;
};

float center_freq;
float phase;
complex float baseband[SCF_SYM_LEN * 2];
complex float fir_tail[SCF_FIR_LEN_RF];

static complex float zadoff_chu_template[FFT_LEN];
static complex float cfo_vector[CFO_STEPS][SCF_SYM_LEN];
static struct rx_chain rx_chain[RX_CHAINS];
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
    for (size_t i = 0; i < FFT_LEN; i++) {
        if (i < SCF_SYM_LEN) {
            fft_buf[i] = in_td[i];
        } else {
            fft_buf[i] = 0.0f;
        }
    }
    fftwf_execute(fft_plan);
    for (size_t i = 0; i < FFT_LEN; i++) {
        fft_buf[i] = fft_buf[i] * template_fd[i];
    }
    fftwf_execute(ifft_plan);
    for (size_t i = 0; i < FFT_LEN; i++) {
        out_td[i] = fft_buf[i] / (float) FFT_LEN;
    }
}

static void demodulate(struct rx_chain *c)
{
    complex float source_cfo[SCF_SYM_LEN] = {0};
    complex float correlated_td[FFT_LEN] = {0};

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        source_cfo[i] = c->source[i] * cfo_vector[c->cfo_index][i];
    }

    correlate_against_template(correlated_td, source_cfo, zadoff_chu_template);

    float max_peak_value = 0.0f;
    size_t max_peak_index = 0;
    for (size_t l = 0; l < SCF_SYM_LEN; l++) {
        int peak1_idx = l;
        int peak2_idx = l - SCF_SYM_LEN + FFT_LEN;
        complex float peak_value_complex = correlated_td[peak1_idx] + correlated_td[peak2_idx];
        float peak_value = peak_value_complex * conjf(peak_value_complex);

        if (peak_value > max_peak_value) {
            max_peak_value = peak_value;
            max_peak_index = l;
        }
    }

    float max_symbol_peak_value = 0.0f;
    size_t max_symbol_index = 0;
    for (uint32_t symbol = 0; symbol < 256; symbol++) {
        size_t symbol_offset = symbol * SCF_DEC_RATIO;
        size_t new_lag = (c->prev_corr_peak_index - symbol_offset + SCF_SYM_LEN) % SCF_SYM_LEN;

        int peak1_idx = new_lag;
        int peak2_idx = new_lag - SCF_SYM_LEN + FFT_LEN;

        complex float peak_value_complex = correlated_td[peak1_idx] + correlated_td[peak2_idx];
        float peak_value = conjf(peak_value_complex) * peak_value_complex;

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
    complex float filtered[SCF_SYM_LEN] = {0};

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        upconverted[i * SCF_DEC_RATIO] = scf_packet_zadoff_chu_sequence[i];
    }

    scf_filter_rf(filtered, upconverted, fir_tail);
    for (size_t i = 0; i < SCF_FIR_LEN_RF; i++) {
        filtered[i] += fir_tail[i];
    }
    for (size_t i = 0; i < FFT_LEN; i++) {
        if (i < SCF_SYM_LEN) {
            fft_buf[i] = filtered[i];
        } else {
            fft_buf[i] = 0.0f;
        }
    }
    fftwf_execute(fft_plan);
    for (size_t i = 0; i < FFT_LEN; i++) {
        zadoff_chu_template[i] = conjf(fft_buf[i]);
    }
}

static void generate_cfo_vector(void)
{
    for (size_t i = 0; i < CFO_STEPS; i++) {
        float cfo = ((float) SCF_BB_SRATE / SCF_BB_SYM_LEN) * (float) i / (float) CFO_STEPS;
        float phase_step = 2.0f * M_PI * cfo / (float) SCF_SRATE;
        float phase = 0.0f;
        for (size_t j = 0; j < SCF_SYM_LEN; j++) {
            phase += phase_step;
            while (phase > 2.0f * M_PI) {
                phase -= 2.0f * M_PI;
            }
            cfo_vector[i][j] = cexpf(-I * phase);
        }
    }
}

void scf_rx_init(float freq)
{
    symbol_counter = 0;
    symbol_skip = 0;

    for (size_t ci = 0; ci < CFO_STEPS; ci++) {
        for (size_t spi = 0; spi < SYM_PHASES; spi++) {
            size_t index = ci * SYM_PHASES + spi;
            rx_chain[index].source = &baseband[spi * SCF_SYM_LEN / SYM_PHASES];
            rx_chain[index].cfo_index = ci;

            // To prevent Reed Solomon decoder from decoding the symbol prematurely
            getrandom(rx_chain[index].decoded_symbol, sizeof(rx_chain[index].decoded_symbol), 0);
        }
    }

    center_freq = freq;
    phase = 0.0f;

    if (!initialized) {
        fft_buf = fftwf_alloc_complex(FFT_LEN);
        assert(fft_buf);

        fft_plan = fftwf_plan_dft_1d(
            FFT_LEN,
            fft_buf,
            fft_buf,
            FFTW_FORWARD,
            FFTW_ESTIMATE
        );
        ifft_plan = fftwf_plan_dft_1d(
            FFT_LEN,
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
    generate_cfo_vector();
}

void scf_rx(scf_rx_result *result, float signal[SCF_SYM_LEN])
{
    memset(result, 0, sizeof(scf_rx_result));
    result->symbol_counter = symbol_counter;

    downconvert(signal);

    for (size_t rx_chain_index = 0; rx_chain_index < RX_CHAINS; rx_chain_index++) {
        struct rx_chain *c = &rx_chain[rx_chain_index];
        demodulate(c);

        if (symbol_skip) {
            continue;
        }

        for (size_t pt = 0; pt < SCF_PACKET_TYPES; pt++) {
            size_t packet_len = scf_packet_fec_len[pt] + 1;
            uint8_t packet[packet_len];
            for (size_t i = 0; i < packet_len; i++) {
                size_t pos = (i + decoded_symbol_index - packet_len + 1 + SCF_MAX_FEC_LEN) % SCF_MAX_FEC_LEN;
                packet[i] = c->decoded_symbol[pos];
            }

            if (!result->got_msg && scf_packet_decode(result, packet, pt)) {
                printf("got cfo index %zu\n", c->cfo_index);
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
