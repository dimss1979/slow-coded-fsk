#ifndef __SCF_H
#define __SCF_H

#include <complex.h>
#include <stddef.h>
#include <stdint.h>

#define SCF_TONES 300
#define SCF_BB_SRATE 2000 /* samples per second */
#define SCF_BB_SYM_LEN 400 /* samples */
#define SCF_DEC_RATIO 4
#define SCF_SRATE (SCF_BB_SRATE * SCF_DEC_RATIO)
#define SCF_SYM_LEN (SCF_BB_SYM_LEN * SCF_DEC_RATIO)
#define SCF_FEC_LEN 8 /* symbols */
#define SCF_PACKET_TYPES 2 /* How many different lengths */
#define SCF_MSG_FEC_MAX 255 /* bytes */

#define SCF_PKT_MAX (SCF_SYNC_LEN + SCF_FEC_LEN * SCF_MSG_FEC_MAX) /* symbols */
#define SCF_SYNC_LEN (5 * SCF_FEC_LEN)
#define SCF_FIR_LEN_RF 369 /* taps */

/* Public API */

void scf_packet_init(uint64_t seed);
size_t scf_encode_packet(uint32_t *packet, uint8_t *msg, size_t msg_len);
void scf_tx_init(float freq);
void scf_tx(float *passband, uint32_t symbol, float gain);
void scf_rx_init(float freq);
size_t scf_rx_symbol(uint8_t *msg, float *signal);

/* Internal */

extern size_t scf_packet_raw_len[SCF_PACKET_TYPES];
extern size_t scf_packet_fec_len[SCF_PACKET_TYPES];
extern void *scf_packet_rs_code[SCF_PACKET_TYPES];
extern uint32_t scf_packet_sync_vector[SCF_PACKET_TYPES][SCF_SYNC_LEN];
extern uint32_t scf_data_code[256][SCF_FEC_LEN];
extern uint8_t scf_data_scrambler[SCF_MSG_FEC_MAX];

void scf_filter_init(void);
void scf_filter_rf(float _Complex *out, float _Complex *in, float _Complex *tail);

#endif
