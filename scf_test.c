#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
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
#include <stdbool.h>

#include "scf.h"

#define MSG_RAW_LEN 107 /* bytes */
#define MSG_FEC_LEN 129 /* bytes */
#define START_OFFSET 20 /* symbols */
#define CARRIER_FREQ 1900.0f /* Hz*/

#define PKT_LEN (SCF_SYNC_LEN + SCF_FEC_LEN * MSG_FEC_LEN) /* symbols */
#define RX_SIGNAL_LEN ((START_OFFSET + PKT_LEN + START_OFFSET) * SCF_SYM_LEN) /* samples */
#define TX_SIGNAL_LEN (PKT_LEN * SCF_SYM_LEN) /* samples */

gsl_rng *rng;

uint8_t tx_message[MSG_RAW_LEN];
uint32_t tx_packet[SCF_PKT_MAX];
float *tx_signal;
float *noise;
float *rx_signal;
uint8_t rx_message[MSG_RAW_LEN];

int wrong_len = 0;
int wrong_msg = 0;

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
        ssize_t rv = write(f, &signal[i], sizeof(signal[i]));
        if (rv != sizeof(signal[i])) {
            break;
        }
    }

    close(f);
}

void generate_message(uint8_t *msg)
{
    for (size_t i = 0; i < MSG_RAW_LEN; i++) {
        msg[i] = gsl_rng_get(rng);
    }
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

bool run_test(void)
{
    memset(tx_signal, 0, RX_SIGNAL_LEN * sizeof(tx_signal[0]));
    memset(noise, 0, RX_SIGNAL_LEN * sizeof(noise[0]));

    scf_packet_init(gsl_rng_get(rng));

    // TX

    size_t tx_rand_delay = gsl_rng_uniform_int(rng, SCF_SYM_LEN) + START_OFFSET * SCF_SYM_LEN;
    float tx_rand_freq_offset = (gsl_rng_uniform(rng) - 0.5f) * 2.0f * 200.0f;

    scf_tx_init(CARRIER_FREQ + tx_rand_freq_offset, (float) SCF_DEC_RATIO);
    generate_message(tx_message);
    bool tx_started = scf_tx_start(tx_message, MSG_RAW_LEN);
    assert(tx_started);

    size_t signal_i = tx_rand_delay;
    for (size_t i = 0; i < PKT_LEN; i++) {
        bool tx_ok = scf_tx(&tx_signal[signal_i]);
        assert(tx_ok);
        signal_i += SCF_SYM_LEN;
    }
    bool tx_ok = scf_tx(&tx_signal[0]);
    assert(!tx_ok);

    //dump_signal("dump_tx_signal.raw", tx_signal, RX_SIGNAL_LEN);
    //exit(1);

    float signal_power = measure_power(&tx_signal[tx_rand_delay], TX_SIGNAL_LEN);

    // Channel

    add_awgn(noise, RX_SIGNAL_LEN, 10.8f);
    add_cw_interferer(noise, RX_SIGNAL_LEN, CARRIER_FREQ, 0.0001f);
    float noise_power = measure_power(&noise[tx_rand_delay], TX_SIGNAL_LEN);

    float snr = 10.0f * log10f(signal_power / noise_power);
    const float user_bitrate = (8.0f * MSG_RAW_LEN * SCF_SRATE) / ((float)PKT_LEN * SCF_SYM_LEN);
    float channel_bandwidth = (float) SCF_SRATE / 2.0f;
    float ebno = snr - 10.0f * log10f(user_bitrate / channel_bandwidth);
    printf(" snr: %f dB\n", snr);
    printf("ebno: %f dB\n", ebno);

    add_signal(rx_signal, tx_signal, noise, RX_SIGNAL_LEN);
    //dump_signal("dump_rx_signal.raw", rx_signal, SIGNAL_LEN);

    // RX

    scf_rx_init(CARRIER_FREQ);
    bool message_is_decoded = false;
    uint64_t decoder_time_start = get_msec();

    for (
        size_t i = 0;
        i < RX_SIGNAL_LEN - SCF_SYM_LEN;
        i += SCF_SYM_LEN
    ) {
        scf_rx_result result;
        scf_rx(&result, &rx_signal[i]);

        if (result.got_sync) {
            printf(" +++ Sync at %zu weight %f cfo %f phase %zu, packet type %zu\n",
                result.symbol_counter, result.sync_weight, result.sync_cfo, result.sync_phase,
                result.sync_packet_type
            );
        }

        if (result.got_msg) {
            if (result.msg_len != MSG_RAW_LEN) {
                printf(" !!! Wrong message length %zu\n", result.msg_len);
                wrong_len++;
                continue;
            }
            if (memcmp(result.msg, tx_message, MSG_RAW_LEN)) {
                printf(" !!! Wrong message\n");
                wrong_msg++;
                continue;
            }
            message_is_decoded = true;
        }
    }

    uint64_t decoder_time = get_msec() - decoder_time_start;
    printf("%s received in %" PRIu64 " msec\n\n", message_is_decoded ? "    " : " NOT", decoder_time);

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
    printf("Wrong message length: %i\n", wrong_len);
    printf("Wrong message: %i\n", wrong_msg);
}
