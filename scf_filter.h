#ifndef __SCF_FILTER_H
#define __SCF_FILTER_H

#include <complex.h>

#define SCF_FIR_LEN_RF 461 /* taps */

void scf_filter_init(void);
void scf_filter_rf(complex float *out, complex float *in, complex float *tail);

#endif
