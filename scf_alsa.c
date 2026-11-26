/* Bridge between SCF audio packets on ALSA and asft-compatible COBS frames
 * on a pseudo-terminal.
 *
 * The bridge is packet length agnostic: it forwards every decoded frame and
 * every decoded SCF message as-is, without validating their length or
 * contents. The upper layer owns packet validation.
 */

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <poll.h>
#include <pthread.h>
#include <pty.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>
#include <alsa/asoundlib.h>

#include "scf.h"
#include "scf_cobs.h"

#define DEFAULT_DEVICE "default"
#define DEFAULT_SYMLINK "/tmp/scf_alsa"
#define DEFAULT_FREQ 1900.0f
#define DEFAULT_SEED 777ULL

#define ALSA_PERIOD SCF_SYM_LEN
#define ALSA_BUFFERS 4
#define PACKET_QUEUE_LEN 4
#define PTY_READ_BUF 256

typedef struct {
    uint8_t data[SCF_MSG_FEC_MAX];
    size_t len;
} packet_t;

/* Only accessed from pty_thread(), so no locking is needed. */
typedef struct {
    packet_t packets[PACKET_QUEUE_LEN];
    size_t head;
    size_t tail;
    size_t count;
} packet_queue_t;

typedef struct {
    snd_pcm_t *capture;
    snd_pcm_t *playback;
    packet_queue_t queue;
} app_t;

static app_t app;
static volatile sig_atomic_t stop_requested;

static int pty_master_fd = -1;
static int pty_slave_fd = -1;
static const char *symlink_path;

static void on_signal(int sig)
{
    (void)sig;
    stop_requested = 1;
}

static void install_signal_handlers(void)
{
    struct sigaction sa;

    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = on_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0; /* do not restart interrupted ALSA/PTY waits */
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);
}

static void usage(FILE *out, const char *prog)
{
    fprintf(out, "Usage: %s [options]\n", prog);
    fprintf(out, "  -d DEVICE   set both capture and playback device (default: %s)\n",
            DEFAULT_DEVICE);
    fprintf(out, "  -C DEVICE   capture device\n");
    fprintf(out, "  -P DEVICE   playback device\n");
    fprintf(out, "  -f HZ       SCF carrier frequency (default: %.0f)\n",
            (double)DEFAULT_FREQ);
    fprintf(out, "  -s SEED     SCF code seed, must match the other end (default: %llu)\n",
            DEFAULT_SEED);
    fprintf(out, "  -l PATH     PTY symlink path (default: %s)\n", DEFAULT_SYMLINK);
    fprintf(out, "  -h          show this help\n");
}

static void packet_queue_init(packet_queue_t *q)
{
    memset(q, 0, sizeof(*q));
}

static bool packet_queue_push(packet_queue_t *q, const uint8_t *data, size_t len)
{
    packet_t *p;

    if (q->count == PACKET_QUEUE_LEN)
        return false;

    p = &q->packets[q->tail];
    memcpy(p->data, data, len);
    p->len = len;
    q->tail = (q->tail + 1) % PACKET_QUEUE_LEN;
    q->count++;

    return true;
}

static bool packet_queue_pop(packet_queue_t *q, packet_t *p)
{
    if (q->count == 0)
        return false;

    *p = q->packets[q->head];
    q->head = (q->head + 1) % PACKET_QUEUE_LEN;
    q->count--;

    return true;
}

static int pcm_open(snd_pcm_t **pcm, const char *name, snd_pcm_stream_t stream)
{
    snd_pcm_hw_params_t *params;
    snd_pcm_uframes_t period;
    snd_pcm_uframes_t buffer;
    unsigned int rate;
    int err;

    err = snd_pcm_open(pcm, name, stream, 0);
    if (err < 0) {
        fprintf(stderr, "Cannot open %s PCM %s: %s\n",
                stream == SND_PCM_STREAM_CAPTURE ? "capture" : "playback",
                name, snd_strerror(err));
        return -1;
    }

    snd_pcm_hw_params_alloca(&params);

    err = snd_pcm_hw_params_any(*pcm, params);
    if (err < 0)
        goto fail;
    err = snd_pcm_hw_params_set_access(*pcm, params, SND_PCM_ACCESS_RW_INTERLEAVED);
    if (err < 0)
        goto fail;
    err = snd_pcm_hw_params_set_format(*pcm, params, SND_PCM_FORMAT_FLOAT);
    if (err < 0)
        goto fail;
    err = snd_pcm_hw_params_set_channels(*pcm, params, 1);
    if (err < 0)
        goto fail;

    rate = SCF_SRATE;
    err = snd_pcm_hw_params_set_rate_near(*pcm, params, &rate, NULL);
    if (err < 0)
        goto fail;
    if (rate != SCF_SRATE) {
        fprintf(stderr, "ALSA device %s does not support %d Hz (got %u Hz)\n",
                name, SCF_SRATE, rate);
        err = -EINVAL;
        goto fail;
    }

    period = ALSA_PERIOD;
    err = snd_pcm_hw_params_set_period_size_near(*pcm, params, &period, NULL);
    if (err < 0)
        goto fail;
    buffer = ALSA_PERIOD * ALSA_BUFFERS;
    err = snd_pcm_hw_params_set_buffer_size_near(*pcm, params, &buffer);
    if (err < 0)
        goto fail;

    err = snd_pcm_hw_params_get_period_size(params, &period, NULL);
    if (err < 0)
        goto fail;
    if (period != ALSA_PERIOD) {
        fprintf(stderr, "ALSA device %s: period %lu instead of %d frames\n",
                name, (unsigned long)period, ALSA_PERIOD);
        err = -EINVAL;
        goto fail;
    }

    err = snd_pcm_hw_params(*pcm, params);
    if (err < 0)
        goto fail;
    err = snd_pcm_prepare(*pcm);
    if (err < 0)
        goto fail;

    fprintf(stderr, "ALSA %s: %s, %u Hz, %lu frames per period\n",
            stream == SND_PCM_STREAM_CAPTURE ? "capture" : "playback",
            name, rate, (unsigned long)period);
    return 0;

fail:
    fprintf(stderr, "ALSA setup failed for %s: %s\n", name, snd_strerror(err));
    snd_pcm_close(*pcm);
    *pcm = NULL;
    return -1;
}

static int pcm_recover(snd_pcm_t *pcm, snd_pcm_sframes_t err)
{
    int ret;

    if (err == -EINTR)
        return 0;
    if (err != -EPIPE && err != -ESTRPIPE)
        return -1;

    ret = snd_pcm_recover(pcm, (int)err, 1);
    if (ret < 0) {
        fprintf(stderr, "ALSA recovery failed: %s\n", snd_strerror(ret));
        return -1;
    }

    return 0;
}

static int pcm_read_full(snd_pcm_t *pcm, float *buf, snd_pcm_uframes_t frames)
{
    snd_pcm_uframes_t done = 0;

    while (done < frames) {
        snd_pcm_sframes_t n = snd_pcm_readi(pcm, buf + done, frames - done);
        if (n > 0) {
            done += (snd_pcm_uframes_t)n;
            continue;
        }
        if (pcm_recover(pcm, n) < 0) {
            fprintf(stderr, "ALSA read error: %s\n", snd_strerror((int)n));
            return -1;
        }
    }

    return 0;
}

static int pcm_write_full(snd_pcm_t *pcm, const float *buf, snd_pcm_uframes_t frames)
{
    snd_pcm_uframes_t done = 0;

    while (done < frames) {
        snd_pcm_sframes_t n = snd_pcm_writei(pcm, buf + done, frames - done);
        if (n > 0) {
            done += (snd_pcm_uframes_t)n;
            continue;
        }
        if (pcm_recover(pcm, n) < 0) {
            fprintf(stderr, "ALSA write error: %s\n", snd_strerror((int)n));
            return -1;
        }
    }

    return 0;
}

static int write_all(int fd, const uint8_t *buf, size_t len)
{
    size_t done = 0;

    while (done < len) {
        ssize_t n = write(fd, buf + done, len - done);
        if (n > 0) {
            done += (size_t)n;
            continue;
        }
        if (n < 0 && errno == EINTR)
            continue;
        if (n < 0 && (errno == EAGAIN || errno == EWOULDBLOCK)) {
            poll(NULL, 0, 10);
            continue;
        }
        perror("PTY write");
        return -1;
    }

    return 0;
}

static void cleanup_pty(void)
{
    if (pty_master_fd >= 0) {
        close(pty_master_fd);
        pty_master_fd = -1;
    }
    if (pty_slave_fd >= 0) {
        close(pty_slave_fd);
        pty_slave_fd = -1;
    }
    if (symlink_path) {
        unlink(symlink_path);
        symlink_path = NULL;
    }
}

static int setup_pty(const char *path)
{
    struct termios t;
    char *slave_name;
    int master;
    int slave;

    if (openpty(&master, &slave, NULL, NULL, NULL) < 0) {
        perror("openpty");
        return -1;
    }

    if (tcgetattr(slave, &t) < 0) {
        perror("tcgetattr");
        goto fail;
    }
    cfmakeraw(&t);
    t.c_cflag |= CLOCAL | CREAD;
    t.c_cc[VMIN] = 1;
    t.c_cc[VTIME] = 0;
    if (tcsetattr(slave, TCSANOW, &t) < 0) {
        perror("tcsetattr");
        goto fail;
    }

    slave_name = ttyname(slave);
    if (!slave_name) {
        fprintf(stderr, "Cannot determine PTY slave name\n");
        goto fail;
    }

    (void)unlink(path);
    if (symlink(slave_name, path) < 0) {
        perror("symlink");
        goto fail;
    }

    pty_master_fd = master;
    pty_slave_fd = slave;
    symlink_path = path;
    fprintf(stderr, "PTY: %s -> %s\n", path, slave_name);
    return 0;

fail:
    close(master);
    close(slave);
    return -1;
}

static void *capture_thread(void *arg)
{
    app_t *a = arg;
    float signal[SCF_SYM_LEN];
    scf_rx_result result;

    while (!stop_requested) {
        if (pcm_read_full(a->capture, signal, SCF_SYM_LEN) < 0)
            break;

        scf_rx(&result, signal);
        if (result.got_msg) {
            uint8_t frame[SCF_COBS_FRAME_MAX];
            size_t frame_len = scf_cobs_frame(result.msg, result.msg_len, frame);

            if (write_all(pty_master_fd, frame, frame_len) < 0)
                break;

            fprintf(stderr, "RX: %zu byte packet, %u FEC corrections\n",
                    result.msg_len, result.outer_fec_errors);
        }
    }

    return NULL;
}

static void *pty_thread(void *arg)
{
    app_t *a = arg;
    uint8_t buf[PTY_READ_BUF];
    uint8_t decoded[SCF_COBS_PKT_MAX];
    size_t decoded_len;
    scf_cobs_rx rx;
    float signal[SCF_SYM_LEN];
    packet_t current;
    bool active = false;

    scf_cobs_rx_reset(&rx);
    memset(signal, 0, sizeof(signal));

    while (!stop_requested) {
        struct pollfd pfd;
        ssize_t n;

        pfd.fd = pty_master_fd;
        pfd.events = POLLIN;
        pfd.revents = 0;

        int ready = poll(&pfd, 1, 0);
        if (ready > 0 && (pfd.revents & POLLIN)) {
            n = read(pty_master_fd, buf, sizeof(buf));
            if (n > 0) {
                for (ssize_t i = 0; i < n; i++) {
                    int r = scf_cobs_rx_push(&rx, buf[i], decoded,
                                             sizeof(decoded), &decoded_len);
                    if (r == 1) {
                        if (packet_queue_push(&a->queue, decoded, decoded_len))
                            fprintf(stderr, "PTY: queued %zu byte packet\n", decoded_len);
                        else
                            fprintf(stderr, "PTY: queue full, dropped %zu byte packet\n",
                                    decoded_len);
                    } else if (r < 0) {
                        fprintf(stderr, "PTY: malformed frame, discarded\n");
                    }
                }
            } else if (n < 0 && errno != EINTR && errno != EAGAIN && errno != EWOULDBLOCK) {
                fprintf(stderr, "PTY read error: %s\n", strerror(errno));
                break;
            }
        } else if (ready < 0 && errno != EINTR) {
            fprintf(stderr, "PTY poll error: %s\n", strerror(errno));
            break;
        }

        if (!active) {
            if (packet_queue_pop(&a->queue, &current)) {
                if (!scf_tx_start(current.data, current.len)) {
                    fprintf(stderr, "TX: cannot start %zu byte packet\n", current.len);
                } else {
                    fprintf(stderr, "TX: sending %zu byte packet\n", current.len);
                    active = true;
                }
            } else {
                memset(signal, 0, sizeof(signal));
            }
        }

        if (active) {
            if (!scf_tx(signal)) {
                active = false;
                memset(signal, 0, sizeof(signal));
            }
        }

        if (pcm_write_full(a->playback, signal, SCF_SYM_LEN) < 0)
            break;
    }

    return NULL;
}

int main(int argc, char **argv)
{
    const char *capture_name = DEFAULT_DEVICE;
    const char *playback_name = DEFAULT_DEVICE;
    const char *symlink = DEFAULT_SYMLINK;
    float freq = DEFAULT_FREQ;
    uint64_t seed = DEFAULT_SEED;
    pthread_t capture_tid;
    pthread_t pty_tid;
    bool capture_started = false;
    bool pty_started = false;
    int opt;

    while ((opt = getopt(argc, argv, "d:C:P:f:s:l:h")) != -1) {
        char *end;

        switch (opt) {
        case 'd':
            capture_name = optarg;
            playback_name = optarg;
            break;
        case 'C':
            capture_name = optarg;
            break;
        case 'P':
            playback_name = optarg;
            break;
        case 'f':
            freq = strtof(optarg, &end);
            if (*end != '\0' || freq <= 0.0f) {
                fprintf(stderr, "Invalid carrier frequency: %s\n", optarg);
                return 1;
            }
            break;
        case 's':
            seed = strtoull(optarg, &end, 0);
            if (*end != '\0') {
                fprintf(stderr, "Invalid seed: %s\n", optarg);
                return 1;
            }
            break;
        case 'l':
            symlink = optarg;
            break;
        case 'h':
            usage(stdout, argv[0]);
            return 0;
        default:
            usage(stderr, argv[0]);
            return 1;
        }
    }

    if (optind != argc) {
        usage(stderr, argv[0]);
        return 1;
    }

    install_signal_handlers();
    signal(SIGPIPE, SIG_IGN);

    scf_packet_init(seed);
    scf_tx_init(freq);
    scf_rx_init(freq);

    if (pcm_open(&app.capture, capture_name, SND_PCM_STREAM_CAPTURE) < 0)
        return 1;
    if (pcm_open(&app.playback, playback_name, SND_PCM_STREAM_PLAYBACK) < 0) {
        snd_pcm_close(app.capture);
        return 1;
    }
    if (setup_pty(symlink) < 0) {
        cleanup_pty();
        snd_pcm_close(app.capture);
        snd_pcm_close(app.playback);
        return 1;
    }

    packet_queue_init(&app.queue);

    if (pthread_create(&capture_tid, NULL, capture_thread, &app) != 0) {
        fprintf(stderr, "Cannot start capture thread\n");
        cleanup_pty();
        snd_pcm_close(app.capture);
        snd_pcm_close(app.playback);
        return 1;
    }
    capture_started = true;

    if (pthread_create(&pty_tid, NULL, pty_thread, &app) != 0) {
        fprintf(stderr, "Cannot start PTY thread\n");
        stop_requested = 1;
        pthread_join(capture_tid, NULL);
        cleanup_pty();
        snd_pcm_close(app.capture);
        snd_pcm_close(app.playback);
        return 1;
    }
    pty_started = true;

    while (!stop_requested)
        poll(NULL, 0, 100);

    if (capture_started)
        pthread_join(capture_tid, NULL);
    if (pty_started)
        pthread_join(pty_tid, NULL);

    cleanup_pty();
    snd_pcm_close(app.capture);
    snd_pcm_close(app.playback);

    return 0;
}
