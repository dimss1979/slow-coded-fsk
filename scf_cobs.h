#ifndef SCF_COBS_H
#define SCF_COBS_H

#include <stddef.h>
#include <stdint.h>

/* Wire format shared with asft: 0x00 | COBS(packet) | 0x00 */
#define SCF_COBS_DELIMITER 0x00
#define SCF_COBS_PKT_MAX 108
#define SCF_COBS_BODY_MAX (SCF_COBS_PKT_MAX + SCF_COBS_PKT_MAX / 254 + 1)
#define SCF_COBS_FRAME_MAX (SCF_COBS_BODY_MAX + 2)

/* COBS encoder/decoder compatible with asft's asft_serial.c. */
size_t scf_cobs_encode(const uint8_t *restrict input, size_t length,
                       uint8_t *restrict output);
size_t scf_cobs_decode(const uint8_t *restrict input, size_t length,
                       uint8_t *restrict output);

/* Build a complete framed packet. Returns the total frame length. */
size_t scf_cobs_frame(const uint8_t *pkt, size_t pkt_len, uint8_t *frame);

typedef struct _scf_cobs_rx {
    uint8_t body[SCF_COBS_BODY_MAX];
    size_t len;
} scf_cobs_rx;

void scf_cobs_rx_reset(scf_cobs_rx *rx);

/* Feed one wire byte.
 * Returns 1 and sets *out_len on a complete, valid frame.
 * Returns 0 while accumulating, -1 on a malformed or oversize frame
 * (the receiver state is reset in that case).
 */
int scf_cobs_rx_push(scf_cobs_rx *rx, uint8_t byte,
                     uint8_t *out, size_t out_cap, size_t *out_len);

#endif
