#include <complex.h>
#include <fftw3.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "scf_private.h"

#define FIR_FFT_LEN 2048
_Static_assert(SCF_SYM_LEN + SCF_FIR_LEN_RF - 1 <= FIR_FFT_LEN, "Invalid FFT length");

/* Low-pass FIR filter for decimation and interpolation.
 * The cutoff is the Nyquist frequency of the decimated baseband.
 */
static fftwf_complex *fft_in_rx, *fft_in_tx, *fft_kernel_fd, *fft_fd, *fft_out;
static fftwf_plan fft, ifft;

static unsigned int initialized;

/* FIR coefficients: Blackman-windowed sinc low pass filter,
 * normalized to unity DC gain.
 */
static void fir_kernel_init(complex float *kernel)
{
    const double bb_srate = (double) SCF_SRATE / (double) SCF_DEC_RATIO;
    const double norm_cutoff = bb_srate / 2.0 / (double) SCF_SRATE;
    const double half = (SCF_FIR_LEN_RF - 1) / 2.0;
    const double window_step = 2.0 * M_PI / (SCF_FIR_LEN_RF - 1);
    double coeff[SCF_FIR_LEN_RF];
    double dc_gain = 0.0;

    for (size_t i = 0; i < SCF_FIR_LEN_RF; i++) {
        double t = (double) i - half;
        double sinc = t == 0.0
            ? 1.0
            : sin(M_PI * 2.0 * norm_cutoff * t) / (M_PI * 2.0 * norm_cutoff * t);
        double window = 0.42
            - 0.5 * cos(window_step * (double) i)
            + 0.08 * cos(2.0 * window_step * (double) i);

        coeff[i] = 2.0 * norm_cutoff * sinc * window;
        dc_gain += coeff[i];
    }

    for (size_t i = 0; i < SCF_FIR_LEN_RF; i++) {
        kernel[i] = coeff[i] / dc_gain;
    }
}

static void scf_filter_internal(complex float *out, complex float *in, fftwf_complex *fft_in)
{
    assert(initialized);

    for (size_t i = 0; i < FIR_FFT_LEN - SCF_SYM_LEN; i++) {
        fft_in[i] = fft_in[i + SCF_SYM_LEN];
    }
    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        fft_in[FIR_FFT_LEN - SCF_SYM_LEN + i] = in[i];
    }

    fftwf_execute_dft(fft, fft_in, fft_fd);

    for (size_t i = 0; i < FIR_FFT_LEN; i++) {
        fft_fd[i] *= fft_kernel_fd[i];
    }

    fftwf_execute(ifft);

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        out[i] = fft_out[FIR_FFT_LEN - SCF_SYM_LEN + i] / (float) FIR_FFT_LEN;
    }
}

void scf_filter_init(void)
{
    if (initialized)
        return;

    fft_in_rx = fftwf_alloc_complex(FIR_FFT_LEN);
    fft_in_tx = fftwf_alloc_complex(FIR_FFT_LEN);
    fft_kernel_fd = fftwf_alloc_complex(FIR_FFT_LEN);
    fft_fd = fftwf_alloc_complex(FIR_FFT_LEN);
    fft_out = fftwf_alloc_complex(FIR_FFT_LEN);
    assert(fft_in_rx && fft_in_tx && fft_kernel_fd && fft_fd && fft_out);

    fft = fftwf_plan_dft_1d(FIR_FFT_LEN, fft_in_rx, fft_fd, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft = fftwf_plan_dft_1d(FIR_FFT_LEN, fft_fd, fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    assert(fft && ifft);

    fir_kernel_init(fft_in_rx);
    for (size_t i = SCF_FIR_LEN_RF; i < FIR_FFT_LEN; i++) {
        fft_in_rx[i] = 0.0f;
    }
    fftwf_execute_dft(fft, fft_in_rx, fft_kernel_fd);

    for (size_t i = 0; i < FIR_FFT_LEN; i++) {
        fft_in_rx[i] = 0.0f;
        fft_in_tx[i] = 0.0f;
    }

    initialized = 1;
}

/* Safe to call before scf_filter_init: the buffers start zeroed there */
void scf_filter_reset_rx(void)
{
    if (initialized) {
        memset(fft_in_rx, 0, FIR_FFT_LEN * sizeof(fft_in_rx[0]));
    }
}

void scf_filter_reset_tx(void)
{
    if (initialized) {
        memset(fft_in_tx, 0, FIR_FFT_LEN * sizeof(fft_in_tx[0]));
    }
}

void scf_filter_rx(complex float *out, complex float *in)
{
    scf_filter_internal(out, in, fft_in_rx);
}

void scf_filter_tx(complex float *out, complex float *in)
{
    scf_filter_internal(out, in, fft_in_tx);
}
