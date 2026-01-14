#ifndef __SCF_RX_H
#define __SCF_RX_H

#include <stdint.h>

void scf_rx_init(float freq);
size_t scf_rx_symbol(uint8_t *msg, float *chip);

#endif
