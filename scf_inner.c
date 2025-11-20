#include <stddef.h>
#include <assert.h>
#include <math.h>
#include <fec.h>

#include "scf.h"
#include "scf_inner.h"

static void *rs_code;

int scf_inner_code[SCF_SYMBOL_M][SCF_CHIPS];

void scf_inner_generate(void)
{
    if (!rs_code) {
        rs_code = init_rs_char(4, 0x13, 1, 1, SCF_CHIPS - 2, 15 - SCF_CHIPS);
    }
    assert(rs_code);

    for (size_t m = 0; m < SCF_SYMBOL_M; m++) {
        uint8_t inner_codeword[SCF_CHIPS];

        inner_codeword[0] = m >> 4;
        inner_codeword[1] = m & 0x0f;

        encode_rs_char(rs_code, &inner_codeword[0], &inner_codeword[2]);

        for (size_t c = 0; c < SCF_CHIPS; c++) {
            scf_inner_code[m][c] = inner_codeword[c];
        }
    }
}
