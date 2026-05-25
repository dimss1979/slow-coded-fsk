#ifndef __SCF_PRIVATE_H
#define __SCF_PRIVATE_H

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

#include "scf.h"

#define SCF_FIR_LEN_RF 369 /* taps */

extern size_t scf_packet_user_len[SCF_PACKET_TYPES];
extern size_t scf_packet_fec_len[SCF_PACKET_TYPES];
extern float _Complex scf_packet_zadoff_chu_sequence[SCF_BB_SYM_LEN];

void scf_filter_init(void);
void scf_filter_rf(float _Complex *out, float _Complex *in, float _Complex *tail);

bool scf_packet_decode(scf_rx_result *result, uint8_t *outer_codeword, size_t outer_codeword_len);

#endif
