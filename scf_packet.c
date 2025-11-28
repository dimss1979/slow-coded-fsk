#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>

#include "scf.h"
#include "scf_packet.h"

static bool initialized = false;

struct scf_packet_type scf_packet_types[SCF_PACKET_TYPES] = {
    {
        .message_len = 10,
        .parity_len = 20,
    },
    {
        .message_len = 100,
        .parity_len = 20,
    }
};

void scf_packet_init(void)
{
    if (initialized)
        return;

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        struct scf_packet_type *pt = &scf_packet_types[i];

        pt->packet_len =  pt->message_len + pt->parity_len;
        assert(pt->packet_len < SCF_SYMBOL_M);
        pt->rs_code = init_rs_char(8, 0x11d, 1, 1, pt->parity_len, 255 - pt->packet_len);
        assert(pt->rs_code);
    }

    initialized = true;
}
