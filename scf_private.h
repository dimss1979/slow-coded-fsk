#ifndef __SCF_PRIVATE_H
#define __SCF_PRIVATE_H

#include <complex.h>
#include <stdbool.h>

#include "scf.h"

#define SCF_FIR_LEN_RF 369 /* taps */

extern size_t scf_packet_raw_len[SCF_PACKET_TYPES];
extern size_t scf_packet_fec_len[SCF_PACKET_TYPES];
extern uint32_t scf_packet_sync_vector[SCF_PACKET_TYPES][SCF_SYNC_LEN];
extern uint32_t scf_inner_code[256][SCF_FEC_LEN];
extern uint8_t scf_data_scrambler[SCF_MSG_FEC_MAX];

void scf_filter_init(void);
void scf_filter_rf(float _Complex *out, float _Complex *in, float _Complex *tail);

bool scf_packet_decode(uint8_t *msg, uint8_t *data_fec_buf, size_t packet_type);

#endif
