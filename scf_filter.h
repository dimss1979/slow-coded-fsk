#ifndef __SCF_FILTER_H
#define __SCF_FILTER_H

#include <complex.h>

#define SCF_FIR_LEN_RF 185 /* taps */

void scf_filter_init(void);
void scf_filter_rf(float _Complex *out, float _Complex *in, float _Complex *tail);

#endif
