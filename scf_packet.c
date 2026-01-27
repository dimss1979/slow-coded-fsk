#include <assert.h>
#include <fec.h>
#include <stddef.h>
#include <stdbool.h>

#include "scf.h"
#include "scf_packet.h"

static bool initialized = false;

size_t scf_packet_raw_len[SCF_PACKET_TYPES] = {
    7,
    107,
};
size_t scf_packet_fec_len[SCF_PACKET_TYPES] = {
    15,
    129,
};
void *scf_packet_rs_code[SCF_PACKET_TYPES];
uint32_t scf_packet_sync_vector[SCF_PACKET_TYPES][SCF_SYNC_LEN];
uint32_t scf_data_code[256][SCF_FEC_LEN];
uint8_t scf_data_scrambler[SCF_MSG_FEC_MAX];

static uint32_t xorshift32(uint32_t *state) {
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;

    return x;
}

static void sync_init(uint32_t *sync_vector, uint32_t seed)
{
    uint32_t xorshift32_state = seed;
    for (size_t i = 0; i < SCF_SYNC_LEN; i++) {
        sync_vector[i] = xorshift32(&xorshift32_state) % SCF_TONES;
    }
}

void scf_packet_init(void)
{
    if (initialized)
        return;

    uint32_t xorshift32_state = 100500;
    for (size_t i = 0; i < 256; i++) {
        for (size_t j = 0; j < SCF_FEC_LEN; j++) {
            scf_data_code[i][j] = xorshift32(&xorshift32_state) % SCF_TONES;
        }
    }

    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        scf_packet_rs_code[i] = init_rs_char(
            8, 0x11d, 1, 1,
            scf_packet_fec_len[i] - scf_packet_raw_len[i],
            255 - scf_packet_fec_len[i]
        );
        assert(scf_packet_rs_code[i]);

        sync_init(&scf_packet_sync_vector[i][0], i + 12345);
    }

    xorshift32_state = 1;
    for (size_t i = 0; i < SCF_MSG_FEC_MAX; i++) {
        scf_data_scrambler[i] = xorshift32(&xorshift32_state);
    }

    initialized = true;
}
