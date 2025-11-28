#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <assert.h>
#include <sys/random.h>

#include "scf.h"
#include "scf_filter.h"
#include "scf_inner.h"
#include "scf_tx.h"
#include "scf_rx.h"
#include "scf_outer.h"
#include "scf_packet.h"

#define MSG_LEN 100 /* bytes */
#define PKT_LEN 120 /* bytes */
#define CRC_LEN 3 /* bytes */

#define START_OFFSET 2
#define RX_SIGNAL_LEN ((START_OFFSET + 1 + SCF_PREAMBLE + PKT_LEN + 6) * SCF_SYMBOL_LEN) /* samples */
#define TX_SIGNAL_LEN ((SCF_PREAMBLE + PKT_LEN) * SCF_SYMBOL_LEN) /* samples */
#define CARRIER_FREQ 1900.0f /* Hz*/

gsl_rng *rng;

uint8_t preamble[SCF_PREAMBLE];
uint8_t tx_message[MSG_LEN];
uint8_t tx_packet[PKT_LEN];
float *tx_signal;
float *noise;
float *rx_signal;
uint8_t rx_message[MSG_LEN];

uint64_t get_msec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);

    return ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
}

void dump_signal(char *filename, float *signal, size_t len)
{
    int f = open(filename, O_CREAT | O_TRUNC | O_RDWR, S_IRUSR | S_IWUSR);
    if (f < 0) {
        return;
    }

    for (size_t i = 0; i < len; i++) {
        int rv = write(f, &signal[i], sizeof(signal[i]));
        assert(rv > 0);
    }

    close(f);
}

static uint32_t crc24(uint8_t *input, size_t len)
{
    uint32_t crc = 0xB704CE;
    while (len--) {
        crc ^= (*input++) << 16;
        for (size_t i = 0; i < 8; i++) {
            crc <<= 1;
            if (crc & 0x1000000) {
                crc &= 0XFFFFFF;
                crc ^= 0x864CFB;
            }
        }
    }
    return crc & 0xFFFFFF;
}

void generate_message(uint8_t *msg)
{
    for (size_t i = 0; i < MSG_LEN; i++) {
        msg[i] = gsl_rng_get(rng) % SCF_SYMBOL_M;
    }

    uint32_t crc = crc24(msg, MSG_LEN - CRC_LEN);

    msg[MSG_LEN - 3] = crc >> 16;
    msg[MSG_LEN - 2] = crc >> 8;
    msg[MSG_LEN - 1] = crc;
}

void add_awgn(float *out, size_t len, float sigma)
{
    for (size_t i = 0; i < len; i++) {
        out[i] += gsl_ran_gaussian(rng, sigma);
    }
}

void add_cw_interferer(float *out, size_t len, float freq, float amplitude)
{
    float phase = 0.0f;

    for (size_t i = 0; i < len; i++) {
        out[i] += amplitude * sinf(phase);
        phase += 2.0f * M_PI * freq * (1.0f / (float) SCF_SRATE);
        if (phase > 2.0f * M_PI) {
            phase -= 2.0f * M_PI;
        }
    }
}

void add_signal(float *out, float *in1, float *in2, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        out[i] = in1[i] + in2[i];
    }
}

float measure_power(float *signal, size_t len)
{
    double energy = 0.0f;

    for (size_t i = 0; i < len; i++) {
        energy += signal[i] * signal[i];
    }

    return energy / (double) len;
}

static bool msg_verifier(uint8_t *msg, size_t msg_len)
{
    if (msg_len < CRC_LEN + 1)
        return false;

    uint32_t crc = crc24(msg, msg_len - CRC_LEN);
    uint32_t received_crc = (msg[msg_len - 3] << 16) | (msg[msg_len - 2] << 8) | msg[msg_len - 1];

    if (crc == received_crc)
        return true;
    else
        return false;
}

bool run_test(void)
{
    memset(noise, 0, RX_SIGNAL_LEN * sizeof(noise[0]));

    // Other user (CDMA test)

    float other_user_level = 0.00001f;

    getrandom(preamble, SCF_PREAMBLE, 0);
    generate_message(tx_message);
    assert(scf_outer_encode(tx_packet, tx_message, MSG_LEN) == PKT_LEN);

    size_t tx_rand_delay = (gsl_rng_uniform(rng) + START_OFFSET) * SCF_SYMBOL_LEN;
    float tx_rand_freq_offset = (gsl_rng_uniform(rng) - 0.5f) * 2.0f * 800.0f;

    scf_tx_init(CARRIER_FREQ + tx_rand_freq_offset);
    size_t signal_i = tx_rand_delay;
    for (size_t i = 0; i < SCF_PREAMBLE; i++) {
        scf_tx(&noise[signal_i], preamble[i], other_user_level);
        signal_i += SCF_SYMBOL_LEN;
    }
    for (size_t i = 0; i < PKT_LEN; i++) {
        scf_tx(&noise[signal_i], tx_packet[i], other_user_level);
        signal_i += SCF_SYMBOL_LEN;
    }

    // TX

    getrandom(preamble, SCF_PREAMBLE, 0);
    generate_message(tx_message);
    assert(scf_outer_encode(tx_packet, tx_message, MSG_LEN) == PKT_LEN);

    tx_rand_delay = (gsl_rng_uniform(rng) + START_OFFSET) * SCF_SYMBOL_LEN;
    tx_rand_freq_offset = (gsl_rng_uniform(rng) - 0.5f) * 2.0f * 800.0f;

    scf_tx_init(CARRIER_FREQ + tx_rand_freq_offset);
    signal_i = tx_rand_delay;
    for (size_t i = 0; i < SCF_PREAMBLE; i++) {
        scf_tx(&noise[signal_i], preamble[i], 1.0f);
        signal_i += SCF_SYMBOL_LEN;
    }
    for (size_t i = 0; i < PKT_LEN; i++) {
        scf_tx(&tx_signal[signal_i], tx_packet[i], 1.0f);
        signal_i += SCF_SYMBOL_LEN;
    }
    //dump_signal("dump_tx_signal.raw", tx_signal, RX_SIGNAL_LEN);
    //exit(1);

    for (size_t i = 0; i < RX_SIGNAL_LEN - 100; i++) {
        // Multipath
        //tx_signal[i] += 1.0f * tx_signal[i + 16];
        //tx_signal[i] *= 0.65f;
    }

    float signal_power = measure_power(&tx_signal[tx_rand_delay], TX_SIGNAL_LEN);

    // Channel

    add_awgn(noise, RX_SIGNAL_LEN, 2.74f);
    add_cw_interferer(noise, RX_SIGNAL_LEN, CARRIER_FREQ, 0.0001f);
    float noise_power = measure_power(&noise[tx_rand_delay], TX_SIGNAL_LEN);

    float snr = 10.0f * log10f(signal_power / noise_power);
    const float user_bitrate = (8.0f * MSG_LEN) / ((SCF_PREAMBLE + PKT_LEN) * SCF_SYMBOL_LEN / SCF_SRATE);
    float channel_bandwidth = (float) SCF_SRATE / 2.0f;
    float ebno = snr - 10.0f * log10f(user_bitrate / channel_bandwidth);
    printf(" snr: %f dB\n", snr);
    printf("ebno: %f dB\n", ebno);

    add_signal(rx_signal, tx_signal, noise, RX_SIGNAL_LEN);
    //dump_signal("dump_rx_signal.raw", rx_signal, SIGNAL_LEN);

    // RX

    scf_rx_init(CARRIER_FREQ, preamble);
    unsigned int symbol_count = 0;
    bool message_is_decoded = false;
    uint64_t decoder_time_start = get_msec();

    for (
        size_t i = 0;
        i < RX_SIGNAL_LEN - SCF_CHIP_LEN + 1;
        i += SCF_CHIP_LEN
    ) {
        size_t received_message_len = scf_rx_chip(rx_message, &rx_signal[i], msg_verifier);
        if (received_message_len) {
            printf("Decoded PSDU %li bytes\n", received_message_len);
            assert(received_message_len == MSG_LEN);
            assert(!memcmp(tx_message, rx_message, MSG_LEN));
            message_is_decoded = true;
            break;
        }
    }

    uint64_t decoder_time = get_msec() - decoder_time_start;
    printf("%s received after %u symbols in %lu msec\n\n", message_is_decoded ? "    " : " NOT", symbol_count, decoder_time);

    return message_is_decoded;
}

int main(int argc, char **argv)
{
    tx_signal = calloc(RX_SIGNAL_LEN, sizeof(tx_signal[0]));
    noise = calloc(RX_SIGNAL_LEN, sizeof(noise[0]));
    rx_signal = calloc(RX_SIGNAL_LEN, sizeof(rx_signal[0]));

    assert(tx_signal);
    assert(noise);
    assert(rx_signal);

    rng = gsl_rng_alloc(gsl_rng_taus2);
    if (!rng) {
        perror("RNG initialization error\n");
        exit(1);
    }
    gsl_rng_set(rng, time(NULL));
    scf_inner_generate();
    scf_packet_init();

    unsigned int decoded_message_cnt = 0;
    const unsigned int test_cnt = 100;
    for (unsigned int i = 0; i < test_cnt; i++) {
        printf("Test %u\n", i);
        if (run_test()) {
            decoded_message_cnt++;
        }
    }
    float reception_probability = (float) decoded_message_cnt / (float) test_cnt;
    printf("Packet reception probability %f\n", reception_probability);
}
