#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <string.h>
#include <fec.h>
#include "scf_private.h"

#define MOD_FILTER_LEN (SCF_BB_SYM_LEN / 20)
#define FREQ_STEP ((float)SCF_BB_SRATE / SCF_BB_SYM_LEN)

static float carrier_freq;
static float carrier_phase;
static float baseband_phase;
static float tx_gain;

static float mod_filter_buf[MOD_FILTER_LEN];
static float mod_filter_kernel[MOD_FILTER_LEN];
static size_t mod_filter_idx;

static uint32_t tx_packet[SCF_PKT_MAX];
static size_t tx_packet_len;
static size_t tx_packet_pos;

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
        int mod_filter_len_i = MOD_FILTER_LEN;
        float mod_filter_len_f = (float) mod_filter_len_i;
        mod_filter_kernel[i] = sinf((M_PI * i) / mod_filter_len_f);
        dc_gain += mod_filter_kernel[i];
    }

    for (size_t i = 0; i < MOD_FILTER_LEN; i++) {
        mod_filter_kernel[i] /= dc_gain;
    }
}

void scf_tx_init(float freq, float gain)
{
    mod_filter_init();
    memset(mod_filter_buf, 0, sizeof(mod_filter_buf));
    mod_filter_idx = 0;
    scf_filter_reset_tx();

    carrier_freq = freq;
    carrier_phase = 0.0f;
    baseband_phase = 0.0f;
    tx_gain = gain;

    tx_packet_len = 0;
    tx_packet_pos = 0;
}

bool scf_tx_start(uint8_t *msg, size_t msg_len)
{
    size_t packet_type = SCF_PACKET_TYPES;
    size_t msg_fec_len = 0;
    void *rs_code = NULL;
    for (size_t i = 0; i < SCF_PACKET_TYPES; i++) {
        if (scf_packet_user_len[i] == msg_len) {
            packet_type = i;
            msg_fec_len = scf_packet_fec_len[i];
            rs_code = scf_packet_rs_code[i];
            break;
        }
    }

    if (packet_type >= SCF_PACKET_TYPES)
        return false;

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
            tx_packet[pos] = scf_inner_code[msg_fec[i]][j];
            pos++;

            sync_cnt--;

            if (!sync_cnt) {
                tx_packet[pos] = scf_packet_sync_vector[packet_type][sync_i];
                pos++;
                sync_i++;
                sync_cnt = block_len;
            }
        }
    }

    tx_packet_len = pos;
    tx_packet_pos = 0;

    return true;
}

bool scf_tx(float signal[SCF_SYM_LEN])
{
    complex float baseband[SCF_SYM_LEN] = {0};
    complex float baseband_filtered[SCF_SYM_LEN] = {0};

    if (tx_packet_pos >= tx_packet_len)
        return false;

    uint32_t symbol = tx_packet[tx_packet_pos];
    tx_packet_pos++;

    float freq = FREQ_STEP * ((float) symbol + 0.5f - (float) SCF_TONES / 2.0f);

    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        baseband[i * SCF_DEC_RATIO] = tx_gain * (complex float) (sinf(baseband_phase) + I * cosf(baseband_phase));

        baseband_phase += 2.0f * M_PI * mod_filter(freq) * (1.0f / SCF_BB_SRATE);
        while (baseband_phase > 2.0f * M_PI) {
            baseband_phase -= 2.0f * M_PI;
        }
        while (baseband_phase < -2.0f * M_PI) {
            baseband_phase += 2.0f * M_PI;
        }
    }

    scf_filter_tx(baseband_filtered, baseband);

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        carrier_phase += 2.0f * M_PI * carrier_freq * (1.0f / SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }

        complex float bb = baseband_filtered[i];
        signal[i] = crealf(bb) * sinf(carrier_phase) - cimagf(bb) * cosf(carrier_phase);
    }

    return true;
}
