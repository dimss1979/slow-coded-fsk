#ifndef __SCF_H
#define __SCF_H

#include <stdint.h>

#define SCF_SYMBOL_M 256 /* fixed for 8 bits per symbol */
#define SCF_CHIPS 8 /* chips per symbol */
#define SCF_SRATE 8000 /* samples per second */
#define SCF_CHIP_LEN 1600 /* samples */
#define SCF_DEC_RATIO 4
#define SCF_PREAMBLE 6 /* bytes */
#define SCF_PREAMBLE_REQUIRED 4 /* bytes */

#define SCF_SYMBOL_LEN (SCF_CHIP_LEN * SCF_CHIPS)
#define SCF_BB_SRATE (SCF_SRATE / SCF_DEC_RATIO)
#define SCF_BB_CHIP_LEN (SCF_CHIP_LEN / SCF_DEC_RATIO)
#define SCF_FREQS (16)
#define CRC_LEN 3

struct scf_soft_symbol {
    int32_t weight;
    uint8_t symbol;
};

#endif
