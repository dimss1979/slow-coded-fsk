#ifndef __SCF_TX_H
#define __SCF_TX_H

#include <complex.h>
#include <stdint.h>

void scf_tx_init(float freq);
size_t scf_encode_packet(uint8_t *packet, uint8_t *msg, size_t msg_raw_len);
void scf_tx(float *passband, uint8_t symbol, float gain);

#endif
