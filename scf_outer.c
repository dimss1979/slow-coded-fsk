#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <endian.h>
#include <string.h>
#include <float.h>

#include "scf.h"
#include "scf_outer.h"
#include "scf_packet.h"

static struct scf_soft_symbol symbol_buf[SCF_SYMBOL_M - 1];

static uint32_t crc24(uint8_t *input, size_t len)
{
    uint32_t crc = 0xB704CE;
    while (len--) {
        crc ^= (*input++) << 16;
        for (size_t i = 0; i < 8; i++) {
            crc <<= 1;
            if (crc & 0x1000000) {
                crc &= 0XFFFFFF;
                crc ^= 0x864CFB;
            }
        }
    }
    return crc & 0xFFFFFF;
}

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

    uint32_t crc = crc24(message, message_len);
    uint8_t *packet_crc = &packet[message_len];
    packet_crc[0] = 0xFF & (crc >> 16);
    packet_crc[1] = 0xFF & (crc >> 8);
    packet_crc[2] = 0xFF & crc;

    uint8_t *packet_parity = &packet[message_len + CRC_LEN];
    encode_rs_char(pt->rs_code, packet, packet_parity);

    return pt->packet_len;
}

size_t scf_outer_decode(uint8_t *message, struct scf_soft_symbol symbol)
{
    memmove(
        &symbol_buf[0],
        &symbol_buf[1],
        sizeof(symbol_buf) - sizeof(symbol_buf[0])
    );
    symbol_buf[SCF_SYMBOL_M - 2] = symbol;

    for (size_t pt_idx = 0; pt_idx < SCF_PACKET_TYPES; pt_idx++) {
        struct scf_packet_type *pt = &scf_packet_types[pt_idx];
        size_t packet_start = SCF_SYMBOL_M - 1 - pt->packet_len;

        int eras_pos[SCF_SYMBOL_M - 1];
        int eras_no_max = pt->parity_len;
        int eras_mark[SCF_SYMBOL_M - 1] = {0};

        for (size_t i = 0; i < eras_no_max; i++) {
            int32_t min_weight = INT32_MAX;
            int min_pos = 0;

            for (size_t j = 0; j < pt->packet_len; j++) {
                if (!eras_mark[j] && symbol_buf[packet_start + j].weight < min_weight) {
                    min_weight = symbol_buf[packet_start + j].weight;
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
                rs_buf[i] = symbol_buf[packet_start + i].symbol;
            }

            int symbol_error_count = decode_rs_char(
                pt->rs_code,
                rs_buf,
                eras_pos_tmp,
                eras_no
            );

            if (symbol_error_count >= 0) {
                uint32_t crc = crc24(rs_buf, pt->message_len);
                uint8_t *packet_crc = &rs_buf[pt->message_len];
                uint32_t crc_received =
                    (packet_crc[0] << 16) | (packet_crc[1] << 8) | packet_crc[2];
                if (crc == crc_received) {
                    memcpy(message, rs_buf, pt->message_len);
                    return pt->message_len;
                }
            }
        }
    }

    return 0;
}

