#ifndef __SCF_PACKET_H
#define __SCF_PACKET_H

#include <stdint.h>

#define SCF_PACKET_TYPES 2

struct scf_packet_type {
    size_t message_len;
    size_t parity_len;
    size_t packet_len;
    void *rs_code;
};

extern struct scf_packet_type scf_packet_types[];

void scf_packet_init(void);

#endif
