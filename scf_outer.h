#ifndef __SCF_OUTER_H
#define __SCF_OUTER_H

#include <stdbool.h>
#include <stdint.h>

#include "scf.h"

size_t scf_outer_encode(uint8_t *packet, uint8_t *message, size_t message_len);
size_t scf_outer_decode(uint8_t *message, struct scf_soft_symbol *symbol_buf, size_t symbol_buf_len, scf_msg_verifier verifier);

#endif
