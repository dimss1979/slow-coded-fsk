#include <stddef.h>
#include <assert.h>
#include <math.h>
#include <fec.h>

#include "scf.h"
#include "scf_inner.h"

static void *rs_code;

int scf_inner_code[SCF_SYMBOL_M][SCF_CHIPS];

static uint32_t xorshift32(uint32_t *state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;

    return x;
}

static uint32_t prng_bounded(uint32_t *state, uint32_t max) {
    return (uint64_t) xorshift32(state) * max >> 32;
}

void permute_uint8(uint8_t *arr, size_t len, uint32_t seed) {
    uint32_t state = seed;

    for (size_t i = len - 1; i > 0; i--) {
        uint32_t j = prng_bounded(&state, i + 1);

        uint8_t tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

void scf_inner_generate(uint32_t seed)
{
    if (!rs_code) {
        rs_code = init_rs_char(4, 0x13, 1, 1, SCF_CHIPS - 2, 15 - SCF_CHIPS);
    }
    assert(rs_code);

    for (size_t m = 0; m < SCF_SYMBOL_M; m++) {
        uint8_t inner_codeword[SCF_CHIPS];

        inner_codeword[0] = m >> 4;
        inner_codeword[1] = m & 0x0f;

        encode_rs_char(rs_code, &inner_codeword[0], &inner_codeword[2]);
        permute_uint8(inner_codeword, SCF_CHIPS, seed);

        for (size_t c = 0; c < SCF_CHIPS; c++) {
            scf_inner_code[m][c] = inner_codeword[c];
        }
    }
}
