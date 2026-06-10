#ifndef __SCF_H
#define __SCF_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* User-configurable parameters
 * Also adjust packet type parameters in scf_packet.c
 */

#define SCF_TONES 300 /* FSK tone count - keep below baseband symbol length */
#define SCF_BB_SRATE 2000 /* Baseband sample rate in samples per second */
#define SCF_BB_SYM_LEN 400 /* Baseband symbol length in samples */
#define SCF_DEC_RATIO 4 /* Decimation ratio */
#define SCF_FEC_LEN 8 /* FSK symbols per inner codeword */
#define SCF_PACKET_TYPES 2 /* How many different lengths */
#define SCF_MSG_FEC_MAX 255 /* Maximum Reed-Solomon outer codeword length */

/* Derived parameters */

#define SCF_SRATE (SCF_BB_SRATE * SCF_DEC_RATIO) /* IO sample rate in samples per second */
#define SCF_SYM_LEN (SCF_BB_SYM_LEN * SCF_DEC_RATIO) /* IO symbol length in samples */
#define SCF_PKT_MAX (SCF_SYNC_LEN + SCF_FEC_LEN * SCF_MSG_FEC_MAX) /* FSK packet length in symbols */
#define SCF_SYNC_LEN (5 * SCF_FEC_LEN) /* Sync vector length in symbols */

typedef struct _scf_rx_result {
    size_t symbol_counter;

    bool got_sync;
    size_t sync_phase;
    size_t sync_packet_type;
    float sync_cfo;
    float sync_weight;

    bool got_msg;
    uint8_t msg[SCF_MSG_FEC_MAX];
    size_t msg_len;
    unsigned int outer_fec_errors;
} scf_rx_result;

/* Always call in the beginning of the program
 * before any other SCF functions. It initializes
 * the packet encoding and decoding tables.
 *
 * Also, call to change parameters:
 *  - seed for the random code and sync vector generation
 *  - carrier frequency
 */
void scf_packet_init(uint64_t seed);
void scf_tx_init(float freq);
void scf_rx_init(float freq);

/* TX: encode a message into a packet of FSK symbol indices */
size_t scf_packet_encode(uint32_t *packet, uint8_t *msg, size_t msg_len);

/* TX: generate a passband FSK signal from an FSK symbol index */
/* Call repeatedly for each FSK symbol in the packet */
void scf_tx(float signal[SCF_SYM_LEN], uint32_t symbol, float gain);

/* RX: decode a message from a passband FSK signal */
/* The result structure is filled according to decoding outcome */
/* Call repeatedly for each passband signal portion of an FSK symbol length */
void scf_rx(scf_rx_result *result, float signal[SCF_SYM_LEN]);

#endif
