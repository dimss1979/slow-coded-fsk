#ifndef __SCF_H
#define __SCF_H

#include <stdint.h>
#include <stdbool.h>

#define SCF_TONES 16
#define SCF_SRATE 8000 /* samples per second */
#define SCF_SYM_LEN 1600 /* samples */
#define SCF_PREAMBLE 40 /* symbols */
#define SCF_PREAMBLE_THR 1400 /* parrots */
#define SCF_HDR_LEN 15 /* symbols */
#define SCF_FEC_LEN 8 /* symbols */
#define SCF_PACKET_TYPES 2 /* How many different lengths */
#define SCF_MSG_FEC_MAX 255 /* bytes */

#define SCF_PKT_MAX (SCF_PREAMBLE + SCF_HDR_LEN + SCF_FEC_LEN * SCF_MSG_FEC_MAX) /* symbols */

typedef bool (*scf_msg_verifier)(uint8_t*, size_t);

#endif
