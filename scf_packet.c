#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "scf_private.h"

static bool rs_code_initialized = false;

/* Packet type parameters.
 * Adjust these to match your needs.
 */

/* User message length 1..253 bytes */
size_t scf_packet_user_len[SCF_PACKET_TYPES] = {
    7,
    107,
};
/* Reed-Solomon outer codeword length 3..255 bytes */
size_t scf_packet_fec_len[SCF_PACKET_TYPES] = {
    15,
    129,
};

static void *scf_packet_rs_code[SCF_PACKET_TYPES];
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

    if (rs_code_initialized)
        return;

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        scf_packet_rs_code[i] = init_rs_char(
            8, 0x11d, 1, 1,
            scf_packet_fec_len[i] - scf_packet_user_len[i],
            255 - scf_packet_fec_len[i]
        );
        assert(scf_packet_rs_code[i]);
    }

    rs_code_initialized = true;
}

size_t scf_packet_encode(uint32_t *packet, uint8_t *msg, size_t msg_len)
{
    size_t packet_type = SCF_PACKET_TYPES;
    size_t msg_fec_len = 0;
    void *rs_code = NULL;
    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        if (scf_packet_user_len[i] == msg_len) {
            packet_type = i;
            msg_fec_len = scf_packet_fec_len[i];
            rs_code = scf_packet_rs_code[i];
            break;
        }
    }
    assert(packet_type < SCF_PACKET_TYPES);

    uint8_t msg_fec[SCF_MSG_FEC_MAX];
    for (size_t i = 0; i < msg_len; i++) {
        msg_fec[i] = msg[i] ^ scf_data_scrambler[i];
    }
    encode_rs_char(rs_code, msg_fec, &msg_fec[msg_len]);

    size_t block_len = (SCF_FEC_LEN * msg_fec_len) / SCF_SYNC_LEN;
    size_t first_block_len = block_len + (SCF_FEC_LEN * msg_fec_len) % SCF_SYNC_LEN;
    size_t pos = 0;
    size_t sync_cnt = first_block_len;
    size_t sync_i = 0;

    for (size_t j = 0; j < SCF_FEC_LEN; j++) {
        for (size_t i = 0; i < msg_fec_len; i++) {
            packet[pos] = scf_inner_code[msg_fec[i]][j];
            pos++;

            sync_cnt--;

            if (!sync_cnt) {
                packet[pos] = scf_packet_sync_vector[packet_type][sync_i];
                pos++;
                sync_i++;
                sync_cnt = block_len;
            }
        }
    }

    return pos;
}

bool scf_packet_decode(scf_rx_result *result, uint8_t *outer_codeword)
{
    size_t packet_type = result->sync_packet_type;
    assert(packet_type < SCF_PACKET_TYPES);

    size_t msg_len = scf_packet_user_len[packet_type];
    int symbol_error_count = decode_rs_char(scf_packet_rs_code[packet_type], outer_codeword, NULL, 0);

    if (symbol_error_count >= 0) {
        for (size_t i = 0; i < msg_len; i++) {
            outer_codeword[i] ^= scf_data_scrambler[i];
        }

        memcpy(result->msg, outer_codeword, msg_len);
        result->msg_len = msg_len;
        result->outer_fec_errors = symbol_error_count;
        result->got_msg = true;

        return true;
    }

    return false;
}
