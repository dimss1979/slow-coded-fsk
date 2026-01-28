#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>

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

void scf_packet_init(uint64_t seed)
{
    uint64_t xorshift64_state = seed + 1000;
    for (size_t i = 0; i < 256; i++) {
        for (size_t j = 0; j < SCF_FEC_LEN; j++) {
            scf_data_code[i][j] = xorshift64(&xorshift64_state) % SCF_TONES;
        }
    }

    xorshift64_state = seed + 2000;
    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        for (size_t j = 0; j < SCF_SYNC_LEN; j++) {
            scf_packet_sync_vector[i][j] = xorshift64(&xorshift64_state) % SCF_TONES;
        }
    }

    xorshift64_state = seed + 3000;
    for (size_t i = 0; i < SCF_MSG_FEC_MAX; i++) {
        scf_data_scrambler[i] = xorshift64(&xorshift64_state);
    }

    if (rs_code_initialized)
        return;

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
