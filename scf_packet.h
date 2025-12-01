#ifndef __SCF_PACKET_H
#define __SCF_PACKET_H

#include <stdint.h>

extern size_t scf_packet_len[];
extern uint8_t scf_header_code[][SCF_HDR_LEN];
extern uint8_t scf_data_code[][SCF_FEC_LEN];

void scf_packet_init(void);

#endif
