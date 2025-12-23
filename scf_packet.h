#ifndef __SCF_PACKET_H
#define __SCF_PACKET_H

#include <stdint.h>

extern size_t scf_packet_raw_len[];
extern size_t scf_packet_fec_len[];
extern void *scf_packet_rs_code[];
extern uint32_t scf_packet_sync_vector[][SCF_PREAMBLE];
extern uint32_t scf_data_code[][SCF_FEC_LEN];
extern uint8_t scf_data_scrambler[];

void scf_packet_init(void);

#endif
