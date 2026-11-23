#ifndef SCF_PRIVATE_H
#define SCF_PRIVATE_H

#include <complex.h>
#include <stdbool.h>

#include "scf.h"

#define SCF_FIR_LEN_RF 369 /* taps */
#define SCF_RS_MAX_LEN 255 /* Longest Reed-Solomon codeword in symbols */

extern size_t scf_packet_user_len[SCF_PACKET_TYPES];
extern size_t scf_packet_fec_len[SCF_PACKET_TYPES];
extern uint32_t scf_packet_sync_vector[SCF_PACKET_TYPES][SCF_SYNC_LEN];
extern uint32_t scf_inner_code[256][SCF_FEC_LEN];
extern uint8_t scf_data_scrambler[SCF_MSG_FEC_MAX];

/* Reed-Solomon outer code: one GF(256) code, hard decision only, no
 * erasures. A codeword is msg_len + parity_len <= SCF_RS_MAX_LEN symbols
 * long and corrects up to parity_len / 2 symbol errors.
 *
 * Appends the parity_len parity symbols to the msg_len message symbols
 * of the codeword. Returns 0, or -1 on invalid parameters.
 */
int scf_rs_encode(uint8_t *codeword, size_t msg_len, size_t parity_len);

/* Corrects the symbol errors of the codeword in place. Returns the
 * number of corrected symbols, 0 if the codeword was already valid, or
 * -1 if it cannot be corrected, in which case it is left unchanged.
 */
int scf_rs_decode(uint8_t *codeword, size_t msg_len, size_t parity_len);

void scf_filter_init(void);
void scf_filter_reset_rx(void);
void scf_filter_reset_tx(void);
void scf_filter_rx(float _Complex *out, float _Complex *in);
void scf_filter_tx(float _Complex *out, float _Complex *in);

#endif
