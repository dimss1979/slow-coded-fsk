#ifndef __SCF_PACKET_H
#define __SCF_PACKET_H

#include <stddef.h>
#include <stdint.h>

#include "scf.h"

extern size_t scf_packet_raw_len[];
extern size_t scf_packet_fec_len[];
extern void *scf_packet_rs_code[];
extern uint32_t scf_packet_sync_vector[][SCF_SYNC_LEN];
extern uint32_t scf_data_code[][SCF_FEC_LEN];
extern uint8_t scf_data_scrambler[];

void scf_packet_init(uint64_t seed);

#endif
