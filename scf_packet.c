#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>

#include "scf.h"
#include "scf_packet.h"

static bool initialized = false;

size_t scf_packet_raw_len[SCF_PACKET_TYPES] = {
    100,
    200,
};
size_t scf_packet_fec_len[SCF_PACKET_TYPES] = {
    120,
    240,
};
void *scf_packet_rs_code[SCF_PACKET_TYPES];
uint8_t scf_header_code[256][SCF_HDR_LEN];
uint8_t scf_data_code[256][SCF_FEC_LEN];


void scf_packet_init(void)
{
    if (initialized)
        return;

    void *rs_header_code = init_rs_char(4, 0x13, 1, 1, SCF_HDR_LEN - 2, 15 - SCF_HDR_LEN);
    void *rs_data_code = init_rs_char(4, 0x13, 1, 1, SCF_FEC_LEN - 2, 15 - SCF_FEC_LEN);

    assert(rs_header_code && rs_data_code);

    for (size_t i = 0; i < 256; i++) {
        scf_header_code[i][0] = (i >> 4) & 0x0F;
        scf_header_code[i][1] = (i >> 0) & 0x0F;
        encode_rs_char(rs_header_code, &scf_header_code[i][0], &scf_header_code[i][2]);

        scf_data_code[i][0] = (i >> 4) & 0x0F;
        scf_data_code[i][1] = (i >> 0) & 0x0F;
        encode_rs_char(rs_data_code, &scf_data_code[i][0], &scf_data_code[i][2]);
    }

    free_rs_char(rs_header_code);
    free_rs_char(rs_data_code);

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        scf_packet_rs_code[i] = init_rs_char(
            8, 0x11d, 1, 1,
            scf_packet_fec_len[i] - scf_packet_raw_len[i],
            255 - scf_packet_fec_len[i]
        );
        assert(scf_packet_rs_code[i]);
    }

    initialized = true;
}
