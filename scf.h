#ifndef __SCF_H
#define __SCF_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/* User-configurable parameters
 * Also adjust packet type parameters in scf_packet.c
 */

#define SCF_BB_SRATE 1000 /* Baseband sample rate in samples per second */
#define SCF_BB_SYM_LEN 997 /* Baseband symbol length in samples */
#define SCF_DEC_RATIO 8 /* Decimation or upsampling ratio */
#define SCF_PACKET_TYPES 2 /* How many different packet lengths */
#define SCF_MAX_FEC_LEN 255 /* Maximum Reed-Solomon outer codeword length */
#define SCF_GUARD_BAND 200.0f /* Guard band in Hz */

/* Derived parameters */

#define SCF_SRATE (SCF_BB_SRATE * SCF_DEC_RATIO) /* IO sample rate in samples per second */
#define SCF_SYM_LEN (SCF_BB_SYM_LEN * SCF_DEC_RATIO) /* IO symbol length in samples */

typedef struct _scf_rx_result {
    size_t symbol_counter;

    bool got_msg;
    uint8_t msg[SCF_MAX_FEC_LEN];
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
void scf_packet_init(void);
void scf_tx_init(float freq);
void scf_rx_init(float freq);

/* TX: encode a message into a packet of ZC symbol indices */
size_t scf_packet_encode(uint32_t *packet, uint8_t *msg, size_t msg_len);

/* TX: generate a passband ZC signal from a ZC symbol index */
/* Call repeatedly for each ZC symbol in the packet */
void scf_tx(float signal[SCF_SYM_LEN], uint32_t symbol, float gain);

/* RX: decode a message from a passband ZC signal */
/* The result structure is filled according to decoding outcome */
/* Call repeatedly for each passband signal portion of a ZC symbol length */
void scf_rx(scf_rx_result *result, float signal[SCF_SYM_LEN]);

#endif
