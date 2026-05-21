#include <math.h>
#include <stddef.h>
#include <assert.h>
#include <string.h>
#include <fec.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_packet.h"
#include "scf_tx.h"

#define MOD_FILTER_LEN (SCF_BB_SYM_LEN / 10)
#define FREQ_STEP ((float)SCF_BB_SRATE / SCF_BB_SYM_LEN)

static float carrier_freq;
static float carrier_phase;
static uint32_t waveform_idx = 0;

static complex float fir_tail[SCF_FIR_LEN_RF];

void scf_tx_init(float freq)
{
    scf_filter_init();

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

    packet[pos] = 0;
    pos++;

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
    waveform_idx = (waveform_idx + 1 + symbol) % SCF_BB_SYM_LEN;

    complex float baseband[SCF_SYM_LEN] = {0};
    for (size_t i = 0; i < SCF_BB_SYM_LEN; i++) {
        baseband[i * SCF_DEC_RATIO] = gain * scf_waveform[waveform_idx][i];
    }

    complex float baseband_filtered[SCF_SYM_LEN] = {0};
    scf_filter_rf(baseband_filtered, baseband, fir_tail);

    for (size_t i = 0; i < SCF_SYM_LEN; i++) {
        carrier_phase += 2.0f * M_PI * carrier_freq * (1.0f / SCF_SRATE);
        while (carrier_phase > 2.0f * M_PI) {
            carrier_phase -= 2.0f * M_PI;
        }

        complex float bb = baseband_filtered[i];
        passband[i] = crealf(bb) * sinf(carrier_phase) - cimagf(bb) * cosf(carrier_phase);
    }
}
