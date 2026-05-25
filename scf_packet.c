#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "scf_private.h"

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

complex float scf_packet_zadoff_chu_sequence[SCF_BB_SYM_LEN];

static bool rs_code_initialized = false;
static void *scf_packet_rs_code[SCF_PACKET_TYPES];

void scf_packet_init(void)
{
    scf_filter_init();

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        scf_packet_zadoff_chu_sequence[i] = cexpf(-I * M_PI * (float) (i * (i + 1)) / (float) SCF_BB_SYM_LEN);
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
    size_t msg_fec_len = 0;
    void *rs_code = NULL;
    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        if (scf_packet_user_len[i] == msg_len) {
            msg_fec_len = scf_packet_fec_len[i];
            rs_code = scf_packet_rs_code[i];
            break;
        }
    }
    assert(msg_fec_len > 0 && rs_code);

    uint8_t msg_fec[SCF_MAX_FEC_LEN];
    for (size_t i = 0; i < msg_len; i++) {
        msg_fec[i] = msg[i];
    }
    encode_rs_char(rs_code, msg_fec, &msg_fec[msg_len]);
    for (size_t i = 0; i < msg_fec_len; i++) {
        packet[i] = msg_fec[i];
    }

    return msg_fec_len;
}

bool scf_packet_decode(scf_rx_result *result, uint8_t *outer_codeword, size_t packet_type)
{
    size_t msg_len = scf_packet_user_len[packet_type];
    int symbol_error_count = decode_rs_char(scf_packet_rs_code[packet_type], outer_codeword, NULL, 0);

    if (symbol_error_count >= 0) {
        memcpy(result->msg, outer_codeword, msg_len);
        result->msg_len = msg_len;
        result->outer_fec_errors = symbol_error_count;
        result->got_msg = true;

        return true;
    }

    return false;
}
