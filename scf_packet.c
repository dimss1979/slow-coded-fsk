#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "scf_private.h"

/* Packet type parameters.
 * Adjust these to match your needs.
 */

/* User message length 1..253 bytes */
size_t scf_packet_user_len[SCF_PACKET_TYPES] = {
    8,
    108,
};
/* Reed-Solomon outer codeword length 3..255 bytes */
size_t scf_packet_fec_len[SCF_PACKET_TYPES] = {
    16,
    130,
};

uint32_t scf_packet_sync_vector[SCF_PACKET_TYPES][SCF_SYNC_LEN];
uint32_t scf_inner_code[256][SCF_FEC_LEN];
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
    scf_filter_init();

    uint64_t xorshift64_state = seed + 1000;
    for (size_t i = 0; i < 256; i++) {
        for (size_t j = 0; j < SCF_FEC_LEN; j++) {
            scf_inner_code[i][j] = xorshift64(&xorshift64_state) % SCF_TONES;
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
}
