/* Test client for scf_alsa's virtual terminal interface.
 *
 * tx: send one 8-byte packet: 4 random bytes plus their little-endian CRC32.
 * rx: continuously receive frames and verify their CRC32.
 */

#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

#include "scf_cobs.h"

#define DEFAULT_PTY "/tmp/scf_alsa"
#define TEST_PKT_LEN 8
#define TEST_RANDOM_LEN 4

static volatile sig_atomic_t stop_requested;
static uint32_t crc32_table[256];

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
    sa.sa_flags = 0; /* do not restart interrupted read()/poll() */
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);
}

static void usage(FILE *out, const char *prog)
{
    fprintf(out, "Usage: %s tx|rx [PTY_PATH]\n", prog);
    fprintf(out, "  tx  send one %d byte packet and exit\n", TEST_PKT_LEN);
    fprintf(out, "  rx  receive packets until interrupted (default path: %s)\n",
            DEFAULT_PTY);
}

static void crc32_init(void)
{
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t c = i;
        for (int k = 0; k < 8; k++)
            c = (c & 1) ? (0xedb88320u ^ (c >> 1)) : (c >> 1);
        crc32_table[i] = c;
    }
}

static uint32_t crc32(const uint8_t *data, size_t len)
{
    uint32_t crc = 0xffffffffu;

    for (size_t i = 0; i < len; i++)
        crc = crc32_table[(crc ^ data[i]) & 0xff] ^ (crc >> 8);

    return crc ^ 0xffffffffu;
}

static bool fill_random(uint8_t *buf, size_t len)
{
    FILE *f = fopen("/dev/urandom", "rb");
    bool ok;

    if (!f)
        return false;
    ok = fread(buf, 1, len, f) == len;
    fclose(f);
    return ok;
}

static void print_hex(const uint8_t *data, size_t len)
{
    for (size_t i = 0; i < len; i++)
        printf("%s%02x", i ? " " : "", data[i]);
}

static int open_pty(const char *path)
{
    struct termios t;
    int fd = open(path, O_RDWR | O_NOCTTY);

    if (fd < 0) {
        perror(path);
        return -1;
    }

    if (tcgetattr(fd, &t) < 0) {
        perror("tcgetattr");
        close(fd);
        return -1;
    }
    cfmakeraw(&t);
    t.c_cflag |= CLOCAL | CREAD;
    t.c_cc[VMIN] = 1;
    t.c_cc[VTIME] = 0;
    if (tcsetattr(fd, TCSANOW, &t) < 0) {
        perror("tcsetattr");
        close(fd);
        return -1;
    }

    return fd;
}

static int run_tx(const char *path)
{
    uint8_t pkt[TEST_PKT_LEN];
    uint8_t frame[SCF_COBS_FRAME_MAX];
    uint32_t crc;
    size_t frame_len;
    size_t done = 0;
    int fd = open_pty(path);

    if (fd < 0)
        return 1;

    if (!fill_random(pkt, TEST_RANDOM_LEN)) {
        srand((unsigned int)time(NULL) ^ (unsigned int)getpid());
        for (size_t i = 0; i < TEST_RANDOM_LEN; i++)
            pkt[i] = (uint8_t)(rand() & 0xff);
    }

    crc = crc32(pkt, TEST_RANDOM_LEN);
    for (size_t i = 0; i < 4; i++)
        pkt[TEST_RANDOM_LEN + i] = (uint8_t)(crc >> (8 * i));

    frame_len = scf_cobs_frame(pkt, sizeof(pkt), frame);

    printf("packet: ");
    print_hex(pkt, sizeof(pkt));
    printf("\npayload: ");
    print_hex(pkt, TEST_RANDOM_LEN);
    printf("\ncrc32: %08x\nframe (%zu bytes): ", crc, frame_len);
    print_hex(frame, frame_len);
    printf("\n");
    fflush(stdout);

    while (done < frame_len) {
        ssize_t n = write(fd, frame + done, frame_len - done);
        if (n > 0) {
            done += (size_t)n;
            continue;
        }
        if (n < 0 && errno == EINTR)
            continue;
        perror("write");
        close(fd);
        return 1;
    }

    close(fd);
    return 0;
}

static int run_rx(const char *path)
{
    uint8_t buf[256];
    uint8_t decoded[SCF_COBS_PKT_MAX];
    size_t decoded_len;
    scf_cobs_rx rx;
    int fd = open_pty(path);

    if (fd < 0)
        return 1;

    scf_cobs_rx_reset(&rx);

    while (!stop_requested) {
        struct pollfd pfd;
        ssize_t n;

        pfd.fd = fd;
        pfd.events = POLLIN;
        pfd.revents = 0;

        int ready = poll(&pfd, 1, 100);
        if (ready < 0) {
            if (errno == EINTR)
                continue;
            perror("poll");
            break;
        }
        if (ready == 0 || !(pfd.revents & POLLIN))
            continue;

        n = read(fd, buf, sizeof(buf));
        if (n < 0) {
            if (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK)
                continue;
            perror("read");
            break;
        }
        if (n == 0)
            continue;

        for (ssize_t i = 0; i < n; i++) {
            time_t now = time(NULL);
            struct tm *tm = localtime(&now);
            char stamp[16];
            int r = scf_cobs_rx_push(&rx, buf[i], decoded,
                                     sizeof(decoded), &decoded_len);

            if (r != 1) {
                if (r < 0)
                    printf("malformed frame, discarded\n");
                continue;
            }

            if (tm)
                strftime(stamp, sizeof(stamp), "%H:%M:%S", tm);
            else
                snprintf(stamp, sizeof(stamp), "--:--:--");

            printf("%s frame (%zu bytes): ", stamp, decoded_len);
            print_hex(decoded, decoded_len);

            if (decoded_len == TEST_PKT_LEN) {
                uint32_t expected = crc32(decoded, TEST_RANDOM_LEN);
                uint32_t got = 0;
                for (size_t j = 0; j < 4; j++)
                    got |= (uint32_t)decoded[TEST_RANDOM_LEN + j] << (8 * j);
                printf("  CRC %s", expected == got ? "OK" : "MISMATCH");
            } else {
                printf("  wrong length (expected %d)", TEST_PKT_LEN);
            }
            printf("\n");
            fflush(stdout);
        }
    }

    close(fd);
    return 0;
}

int main(int argc, char **argv)
{
    const char *path = DEFAULT_PTY;

    if (argc < 2 || argc > 3) {
        usage(stderr, argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "tx") != 0 && strcmp(argv[1], "rx") != 0) {
        usage(stderr, argv[0]);
        return 1;
    }
    if (argc == 3)
        path = argv[2];

    install_signal_handlers();
    crc32_init();

    if (strcmp(argv[1], "tx") == 0)
        return run_tx(path);
    return run_rx(path);
}
