#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <endian.h>
#include <string.h>
#include <float.h>

#include "scf.h"
#include "scf_outer.h"
#include "scf_packet.h"

size_t scf_outer_encode(uint8_t *packet, uint8_t *message, size_t message_len)
{
    struct scf_packet_type *pt = NULL;

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        if (scf_packet_types[i].message_len == message_len) {
            pt = &scf_packet_types[i];
            break;
        }
    }

    if (!pt)
        return 0;

    memcpy(packet, message, message_len);
    encode_rs_char(pt->rs_code, packet, &packet[message_len]);

    return pt->packet_len;
}

size_t scf_outer_decode(uint8_t *message, struct scf_soft_symbol *symbol_buf, size_t symbol_buf_len, scf_msg_verifier verifier)
{
    for (size_t pt_idx = 0; pt_idx < SCF_PACKET_TYPES; pt_idx++) {
        struct scf_packet_type *pt = &scf_packet_types[pt_idx];

        if (symbol_buf_len != pt->packet_len)
            continue;

        int eras_pos[SCF_SYMBOL_M - 1];
        int eras_no_max = pt->parity_len;
        int eras_mark[SCF_SYMBOL_M - 1] = {0};

        for (size_t i = 0; i < eras_no_max; i++) {
            int32_t min_weight = INT32_MAX;
            int min_pos = 0;

            for (size_t j = 0; j < pt->packet_len; j++) {
                if (!eras_mark[j] && symbol_buf[j].weight < min_weight) {
                    min_weight = symbol_buf[j].weight;
                    min_pos = j;
                }
            }

            eras_mark[min_pos] = 1;
            eras_pos[i] = min_pos;
        }

        for (size_t eras_no = 0; eras_no <= eras_no_max; eras_no += 2) {
            int eras_pos_tmp[SCF_SYMBOL_M - 1];
            memcpy(eras_pos_tmp, eras_pos, sizeof(eras_pos_tmp));

            uint8_t rs_buf[SCF_SYMBOL_M - 1];
            for (size_t i = 0; i < pt->packet_len; i++) {
                rs_buf[i] = symbol_buf[i].symbol;
            }

            int symbol_error_count = decode_rs_char(
                pt->rs_code,
                rs_buf,
                eras_pos_tmp,
                eras_no
            );

            if (symbol_error_count >= 0) {
                if (verifier(rs_buf, pt->message_len)) {
                    memcpy(message, rs_buf, pt->message_len);
                    return pt->message_len;
                }
            }
        }
    }

    return 0;
}

