#ifndef __SCF_H
#define __SCF_H

#include <stdint.h>
#include <stdbool.h>

#define SCF_TONES 200
#define SCF_SRATE 8000 /* samples per second */
#define SCF_SYM_LEN 1600 /* samples */
#define SCF_FEC_LEN 8 /* symbols */
#define SCF_PACKET_TYPES 2 /* How many different lengths */
#define SCF_MSG_FEC_MAX 255 /* bytes */

#define SCF_PKT_MAX (SCF_SYNC_LEN + SCF_FEC_LEN * SCF_MSG_FEC_MAX) /* symbols */
#define SCF_SYNC_LEN (5 * SCF_FEC_LEN)

#endif
