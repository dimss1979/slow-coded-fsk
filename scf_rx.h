#ifndef __SCF_RX_H
#define __SCF_RX_H

#include <stdint.h>

void scf_rx_init(float freq, uint8_t *preamble);
size_t scf_rx_chip(uint8_t *msg, float *chip, scf_msg_verifier verifier);

#endif
