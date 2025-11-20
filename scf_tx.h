#ifndef __SCF_TX_H
#define __SCF_TX_H

#include <complex.h>
#include <stdint.h>

void scf_tx_init(float freq);
void scf_tx(float *passband, uint8_t symbol, float gain);

#endif
