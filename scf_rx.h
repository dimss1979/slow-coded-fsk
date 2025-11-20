#ifndef __SCF_RX_H
#define __SCF_RX_H

#include <stdbool.h>
#include <stdint.h>

void scf_rx_init(float freq, uint8_t *preamble);
bool scf_rx_chip(struct scf_soft_symbol *symbol, float *chip);

#endif
