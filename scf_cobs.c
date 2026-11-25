/* COBS framing ported from asft's asft_serial.c to remain wire-compatible.
 *
 * asft is distributed under the ISC License:
 *
 * Copyright (c) 2026 dimss1979 <dimss1979@yandex.ru>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include "scf_cobs.h"

size_t scf_cobs_encode(const uint8_t *restrict input, size_t length,
                       uint8_t *restrict output)
{
    size_t read_index = 0;
    size_t write_index = 1;
    size_t code_index = 0;
    uint8_t code = 1;

    while (read_index < length) {
        if (input[read_index] == 0) {
            output[code_index] = code;
            code = 1;
            code_index = write_index;
            write_index++;
            read_index++;
        } else {
            output[write_index] = input[read_index];
            write_index++;
            read_index++;
            code++;
            if (code == 0xff) {
                output[code_index] = code;
                code = 1;
                code_index = write_index;
                write_index++;
            }
        }
    }

    output[code_index] = code;
    return write_index;
}

size_t scf_cobs_decode(const uint8_t *restrict input, size_t length,
                       uint8_t *restrict output)
{
    size_t read_index = 0;
    size_t write_index = 0;
    uint8_t code;
    uint8_t i;

    while (read_index < length) {
        code = input[read_index];

        if (read_index + code > length && code != 1)
            return 0;

        read_index++;

        for (i = 1; i < code; i++) {
            output[write_index] = input[read_index];
            write_index++;
            read_index++;
        }

        if (code != 0xff && read_index != length) {
            output[write_index] = 0x00;
            write_index++;
        }
    }

    return write_index;
}

size_t scf_cobs_frame(const uint8_t *pkt, size_t pkt_len, uint8_t *frame)
{
    frame[0] = SCF_COBS_DELIMITER;
    size_t body_len = scf_cobs_encode(pkt, pkt_len, &frame[1]);
    frame[1 + body_len] = SCF_COBS_DELIMITER;
    return body_len + 2;
}

void scf_cobs_rx_reset(scf_cobs_rx *rx)
{
    rx->len = 0;
}

int scf_cobs_rx_push(scf_cobs_rx *rx, uint8_t byte,
                     uint8_t *out, size_t out_cap, size_t *out_len)
{
    if (byte == SCF_COBS_DELIMITER) {
        if (rx->len == 0)
            return 0;

        size_t decoded_len = scf_cobs_decode(rx->body, rx->len, out);
        rx->len = 0;

        if (decoded_len == 0 || decoded_len > out_cap)
            return -1;

        *out_len = decoded_len;
        return 1;
    }

    if (rx->len >= SCF_COBS_BODY_MAX) {
        rx->len = 0;
        return -1;
    }

    rx->body[rx->len++] = byte;
    return 0;
}
