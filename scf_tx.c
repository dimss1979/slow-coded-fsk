#include <math.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <fec.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_packet.h"
#include "scf_tx.h"

#define MOD_FILTER_LEN (SCF_SYM_LEN / 10)

static float carrier_freq;
static float carrier_phase;
static const float freq_step = (float) SCF_SRATE / SCF_SYM_LEN;

static float mod_filter_buf[MOD_FILTER_LEN];
static float mod_filter_kernel[MOD_FILTER_LEN];
static size_t mod_filter_idx;

static float mod_filter(float x)
{
    mod_filter_buf[mod_filter_idx] = x;
    mod_filter_idx = (mod_filter_idx + 1) % MOD_FILTER_LEN;

    float y = 0.0f;
    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        size_t pos = (mod_filter_idx + i) % MOD_FILTER_LEN;
        y += mod_filter_buf[pos] * mod_filter_kernel[i];
    }

    return y;
}

static void mod_filter_init(void)
{
    float dc_gain = 0.0f;

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        mod_filter_kernel[i] = sin(M_PI * i / MOD_FILTER_LEN);
        dc_gain += mod_filter_kernel[i];
    }

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        mod_filter_kernel[i] /= dc_gain;
    }
}

void scf_tx_init(float freq)
{
    scf_filter_init();
    mod_filter_init();

    carrier_freq = freq;
}

size_t scf_encode_packet(uint32_t *packet, uint8_t *msg, size_t msg_len)
{
    size_t packet_type = -1;
    size_t msg_fec_len = 0;
    void *rs_code = NULL;
    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        if (scf_packet_raw_len[i] == msg_len) {
            packet_type = i;
            msg_fec_len = scf_packet_fec_len[i];
            rs_code = scf_packet_rs_code[i];
            break;
        }
    }
    assert(packet_type >= 0);

    uint8_t msg_fec[SCF_MSG_FEC_MAX];
    for (size_t i = 0; i < msg_len; i++) {
        msg_fec[i] = msg[i] ^ scf_data_scrambler[i];
    }
    encode_rs_char(rs_code, msg_fec, &msg_fec[msg_len]);

    size_t block_len = (SCF_FEC_LEN * msg_fec_len) / SCF_SYNC_LEN;
    size_t first_block_len = block_len + (SCF_FEC_LEN * msg_fec_len) % SCF_SYNC_LEN;
    size_t pos = 0;
    size_t sync_cnt = first_block_len;
    size_t sync_i = 0;

    for (size_t j = 0; j < SCF_FEC_LEN; j++) {
        for (size_t i = 0; i < msg_fec_len; i++) {
            packet[pos] = scf_data_code[msg_fec[i]][j];
            pos++;

            sync_cnt--;

            if (!sync_cnt) {
                packet[pos] = scf_packet_sync_vector[packet_type][sync_i];
                pos++;
                sync_i++;
                sync_cnt = block_len;
            }
        }
    }

    return pos;
}

void scf_tx(float *passband, uint32_t symbol, float gain)
{
    float freq = carrier_freq + freq_step * symbol;

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        passband[i] = gain * sinf(carrier_phase);

        carrier_phase += 2.0f * M_PI * mod_filter(freq) * (1.0f / SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }
    }
}
