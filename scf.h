#ifndef __SCF_H
#define __SCF_H

#include <stdint.h>
#include <stdbool.h>

#define SCF_BB_SRATE 2000 /* samples per second */
#define SCF_BB_SYM_LEN 400 /* samples */
#define SCF_TONES (SCF_BB_SYM_LEN - 1)
#define SCF_DEC_RATIO 4
#define SCF_SRATE (SCF_BB_SRATE * SCF_DEC_RATIO)
#define SCF_SYM_LEN (SCF_BB_SYM_LEN * SCF_DEC_RATIO)
#define SCF_FEC_LEN 8 /* symbols */
#define SCF_PACKET_TYPES 2 /* How many different lengths */
#define SCF_MSG_FEC_MAX 255 /* bytes */

#define SCF_PKT_MAX (SCF_SYNC_LEN + SCF_FEC_LEN * SCF_MSG_FEC_MAX + 1) /* symbols */
#define SCF_SYNC_LEN (5 * SCF_FEC_LEN)

#endif
