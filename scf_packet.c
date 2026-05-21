#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>

#include "scf.h"
#include "scf_packet.h"

static bool rs_code_initialized = false;

size_t scf_packet_raw_len[SCF_PACKET_TYPES] = {
    7,
    107,
};
size_t scf_packet_fec_len[SCF_PACKET_TYPES] = {
    15,
    129,
};
void *scf_packet_rs_code[SCF_PACKET_TYPES];
uint32_t scf_packet_sync_vector[SCF_PACKET_TYPES][SCF_SYNC_LEN];
uint32_t scf_data_code[256][SCF_FEC_LEN];
uint8_t scf_data_scrambler[SCF_MSG_FEC_MAX];
complex float scf_waveform[SCF_BB_SYM_LEN][SCF_BB_SYM_LEN];

static uint64_t xorshift64(uint64_t *state) {
    uint64_t x = *state;

    if (!x)
        x = 1;

    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;

    return x;
}

static void scf_packet_init_data_code(uint64_t seed)
{
    uint64_t xorshift64_state = seed + 1000;

    for (size_t i = 0; i < 256; i++) {
        for (size_t j = 0; j < SCF_FEC_LEN; j++) {
            scf_data_code[i][j] = xorshift64(&xorshift64_state) % SCF_TONES;
        }
    }
}

static void scf_packet_init_sync_vector(uint64_t seed)
{
    uint64_t xorshift64_state = seed + 2000;

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        for (size_t j = 0; j < SCF_SYNC_LEN; j++) {
            scf_packet_sync_vector[i][j] = xorshift64(&xorshift64_state) % SCF_TONES;
        }
    }
}

static void scf_packet_init_data_scrambler(uint64_t seed)
{
    uint64_t xorshift64_state = seed + 3000;

    for (size_t i = 0; i < SCF_MSG_FEC_MAX; i++) {
        scf_data_scrambler[i] = xorshift64(&xorshift64_state);
    }
}

static void scf_packet_init_rs_code(void)
{
    if (rs_code_initialized) {
        return;
    }

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        scf_packet_rs_code[i] = init_rs_char(
            8, 0x11d, 1, 1,
            scf_packet_fec_len[i] - scf_packet_raw_len[i],
            255 - scf_packet_fec_len[i]
        );
        assert(scf_packet_rs_code[i]);
    }

    rs_code_initialized = true;
}

void scf_waveform_init(uint64_t seed)
{
    fftwf_complex *fft_buf = fftwf_alloc_complex(SCF_BB_SYM_LEN);
    assert(fft_buf);
    fftwf_plan fft_plan = fftwf_plan_dft_1d(
        SCF_BB_SYM_LEN,
        fft_buf,
        fft_buf,
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );
    assert(fft_plan);

    uint64_t xorshift64_state = seed + 228264;
    complex float primary_spectrum[SCF_BB_SYM_LEN];

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        uint32_t random_number = xorshift64(&xorshift64_state);
        float magnitude = (float) (random_number % 2);
        float phase = (float) (random_number % 360) * M_PI / 180.0f;

        primary_spectrum[i] = magnitude * (complex float) (sinf(phase) + I * cosf(phase));
    }

    for (size_t w = 0; w < SCF_BB_SYM_LEN; w++) {
        for (size_t s = 0; s < SCF_BB_SYM_LEN; s++) {
            fft_buf[s] = primary_spectrum[(s + w) % SCF_BB_SYM_LEN];
        }
        fft_buf[0] = 0.0f;

        fftwf_execute(fft_plan);

        float energy = 0.0f;
        for (size_t s = 0; s < SCF_BB_SYM_LEN; s++) {
            energy += fft_buf[s] * conjf(fft_buf[s]);
        }
        float rms = sqrtf(energy / SCF_BB_SYM_LEN);

        for (size_t s = 0; s < SCF_BB_SYM_LEN; s++) {
            scf_waveform[w][s] = fft_buf[s] / rms;
        }
    }
}

void scf_packet_init(uint64_t seed)
{
    scf_packet_init_data_code(seed);
    scf_packet_init_sync_vector(seed);
    scf_packet_init_data_scrambler(seed);
    scf_packet_init_rs_code();
    scf_waveform_init(seed);
}
