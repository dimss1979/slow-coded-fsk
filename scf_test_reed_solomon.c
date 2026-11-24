/* Tests for the Reed-Solomon codec in scf_reed_solomon.c
 *
 * Agent: OpenCode, model: MiMo-V2.6-Flash Free
 */

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "scf_private.h"

static unsigned int checks;
static unsigned int failures;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { \
        failures++; \
        printf("FAIL: " __VA_ARGS__); \
        printf("\n"); \
    } \
} while (0)

struct rs_params {
    size_t msg_len;
    size_t parity_len;
};

/* A selection of message and parity lengths, covering the packet types
 * used by SCF and both ends of the supported range */
static const struct rs_params params[] = {
    {1, 2},    /* Smallest correctable code */
    {8, 8},    /* SCF short packet */
    {20, 12},
    {50, 3},   /* Odd number of parity symbols */
    {108, 22}, /* SCF long packet */
    {100, 100},
    {128, 127},
    {253, 2},  /* Longest message */
    {1, 254},  /* Largest parity */
};

#define PARAMS_LEN (sizeof(params) / sizeof(params[0]))

static uint64_t rng_state = 0x2545f4914f6cdd1dULL;

static uint64_t rng_next(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;

    return rng_state;
}

static size_t rng_below(size_t n)
{
    return (size_t) (rng_next() % n);
}

static void fill_random(uint8_t *buf, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        buf[i] = (uint8_t) (rng_next() >> 32);
    }
}

/* Symbol change that never leaves the symbol untouched */
static uint8_t nonzero_delta(void)
{
    uint8_t delta = (uint8_t) (rng_next() >> 32);

    return delta ? delta : 1;
}

static void shuffle(size_t *positions, size_t len)
{
    for (size_t i = len - 1; i > 0; i--) {
        size_t j = rng_below(i + 1);
        size_t tmp = positions[i];
        positions[i] = positions[j];
        positions[j] = tmp;
    }
}

static void make_codeword(uint8_t *codeword, const struct rs_params *p)
{
    fill_random(codeword, p->msg_len);
    scf_rs_encode(codeword, p->msg_len, p->parity_len);
}

/* Independent GF(256) multiply with the primitive polynomial 0x11d */
static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    uint8_t product = 0;

    for (int i = 0; i < 8; i++) {
        if (b & 1)
            product ^= a;

        b >>= 1;

        if (a & 0x80) {
            a <<= 1;
            a ^= 0x1d;
        } else {
            a <<= 1;
        }
    }

    return product;
}

/* Is every syndrome of the codeword zero? */
static bool is_codeword(const uint8_t *codeword, size_t msg_len, size_t parity_len)
{
    size_t len = msg_len + parity_len;
    uint8_t root = 2; /* alpha^1, the first root of the generator */

    for (size_t i = 0; i < parity_len; i++) {
        uint8_t syndrome = 0;

        for (size_t j = 0; j < len; j++)
            syndrome = gf_mul(syndrome, root) ^ codeword[j];

        if (syndrome)
            return false;

        root = gf_mul(root, 2);
    }

    return true;
}

static void test_parameters(void)
{
    uint8_t codeword[SCF_RS_MAX_LEN];

    memset(codeword, 0, sizeof(codeword));

    CHECK(scf_rs_encode(NULL, 4, 4) == -1, "encode without a codeword");
    CHECK(scf_rs_encode(codeword, 0, 4) == -1, "encode without a message");
    CHECK(scf_rs_encode(codeword, 4, 0) == -1, "encode without a parity");
    CHECK(scf_rs_encode(codeword, SCF_RS_MAX_LEN, 1) == -1,
        "encode beyond the maximum codeword length");
    CHECK(scf_rs_encode(codeword, 1, SCF_RS_MAX_LEN) == -1,
        "encode beyond the maximum codeword length");
    CHECK(scf_rs_encode(codeword, SCF_RS_MAX_LEN - 1, 1) == 0,
        "encode of a maximum length codeword");

    CHECK(scf_rs_decode(NULL, 4, 4) == -1, "decode without a codeword");
    CHECK(scf_rs_decode(codeword, 0, 4) == -1, "decode without a message");
    CHECK(scf_rs_decode(codeword, 4, 0) == -1, "decode without a parity");
    CHECK(scf_rs_decode(codeword, SCF_RS_MAX_LEN, 1) == -1,
        "decode beyond the maximum codeword length");
}

static void test_known_parity(void)
{
    /* Parity of an incrementing message as generated with libfec's
     * init_rs_char(8, 0x11d, 1, 1, nroots, pad). Packets that are
     * already on the air must keep decoding, so this must not change.
     */
    static const struct {
        size_t msg_len;
        size_t parity_len;
        uint8_t parity[22];
    } known[] = {
        {7, 8, {
            0x7f, 0xfb, 0x09, 0xf7, 0x78, 0x08, 0x81, 0x07,
        }},
        {107, 22, {
            0x0b, 0x32, 0x8c, 0xf9, 0x34, 0xae, 0xe7, 0x51, 0xf1, 0xb7,
            0xec, 0xec, 0x60, 0x9c, 0x35, 0x02, 0xa2, 0x7b, 0x65, 0xa3,
            0xcf, 0x7a,
        }},
    };

    uint8_t codeword[SCF_RS_MAX_LEN];

    for (size_t i = 0; i < sizeof(known) / sizeof(known[0]); i++) {
        for (size_t j = 0; j < known[i].msg_len; j++) {
            codeword[j] = (uint8_t) j;
        }

        CHECK(scf_rs_encode(codeword, known[i].msg_len, known[i].parity_len) == 0,
            "encode of the %zu symbol message", known[i].msg_len);
        CHECK(!memcmp(&codeword[known[i].msg_len], known[i].parity, known[i].parity_len),
            "parity of the %zu symbol message matches the known one", known[i].msg_len);
    }
}

static void test_roundtrip(void)
{
    uint8_t message[SCF_RS_MAX_LEN];
    uint8_t codeword[SCF_RS_MAX_LEN];
    uint8_t reference[SCF_RS_MAX_LEN];

    for (size_t i = 0; i < PARAMS_LEN; i++) {
        const struct rs_params *p = &params[i];
        size_t len = p->msg_len + p->parity_len;

        fill_random(message, p->msg_len);
        memcpy(codeword, message, p->msg_len);
        CHECK(scf_rs_encode(codeword, p->msg_len, p->parity_len) == 0,
            "encode of a %zu+%zu codeword", p->msg_len, p->parity_len);
        memcpy(reference, codeword, len);

        /* The parity only depends on the message */
        memcpy(codeword, message, p->msg_len);
        scf_rs_encode(codeword, p->msg_len, p->parity_len);
        CHECK(!memcmp(codeword, reference, len),
            "parity of the %zu+%zu codeword is reproducible", p->msg_len, p->parity_len);

        CHECK(scf_rs_decode(codeword, p->msg_len, p->parity_len) == 0,
            "decode of a clean %zu+%zu codeword", p->msg_len, p->parity_len);
        CHECK(!memcmp(codeword, reference, len),
            "clean %zu+%zu codeword stays intact", p->msg_len, p->parity_len);
        CHECK(!memcmp(codeword, message, p->msg_len),
            "message of the clean %zu+%zu codeword stays intact", p->msg_len, p->parity_len);
    }
}

static void test_single_errors(void)
{
    const struct rs_params p = {20, 12};
    uint8_t codeword[SCF_RS_MAX_LEN];
    uint8_t reference[SCF_RS_MAX_LEN];
    size_t len = p.msg_len + p.parity_len;

    make_codeword(codeword, &p);
    memcpy(reference, codeword, len);

    for (size_t i = 0; i < len; i++) {
        uint8_t delta = nonzero_delta();

        memcpy(codeword, reference, len);
        codeword[i] ^= delta;

        int fixed = scf_rs_decode(codeword, p.msg_len, p.parity_len);
        CHECK(fixed == 1, "single error at symbol %zu of %zu corrected, got %d",
            i, len, fixed);
        CHECK(!memcmp(codeword, reference, len),
            "codeword restored after the single error at symbol %zu", i);
    }
}

static void test_correctable_errors(void)
{
    uint8_t codeword[SCF_RS_MAX_LEN];
    uint8_t reference[SCF_RS_MAX_LEN];
    size_t positions[SCF_RS_MAX_LEN];

    for (size_t i = 0; i < PARAMS_LEN; i++) {
        const struct rs_params *p = &params[i];
        size_t len = p->msg_len + p->parity_len;
        size_t correctable = p->parity_len / 2;

        make_codeword(codeword, p);
        memcpy(reference, codeword, len);

        for (size_t j = 0; j < len; j++) {
            positions[j] = j;
        }
        shuffle(positions, len);

        for (size_t errors = 1; errors <= correctable; errors++) {
            memcpy(codeword, reference, len);
            for (size_t k = 0; k < errors; k++) {
                codeword[positions[k]] ^= nonzero_delta();
            }

            int fixed = scf_rs_decode(codeword, p->msg_len, p->parity_len);
            CHECK(fixed == (int) errors,
                "%zu errors of a %zu+%zu codeword corrected, got %d",
                errors, p->msg_len, p->parity_len, fixed);
            CHECK(!memcmp(codeword, reference, len),
                "%zu+%zu codeword restored with %zu errors",
                p->msg_len, p->parity_len, errors);
        }
    }
}

static void test_uncorrectable_errors(void)
{
    uint8_t codeword[SCF_RS_MAX_LEN];
    uint8_t before[SCF_RS_MAX_LEN];
    size_t positions[SCF_RS_MAX_LEN];

    for (size_t i = 0; i < PARAMS_LEN; i++) {
        const struct rs_params *p = &params[i];
        size_t len = p->msg_len + p->parity_len;
        size_t errors = p->parity_len / 2 + 1;

        make_codeword(codeword, p);

        for (size_t j = 0; j < len; j++) {
            positions[j] = j;
        }
        shuffle(positions, len);

        for (size_t j = 0; j < errors; j++) {
            codeword[positions[j]] ^= nonzero_delta();
        }
        memcpy(before, codeword, len);

        int fixed = scf_rs_decode(codeword, p->msg_len, p->parity_len);

        if (fixed == -1) {
            CHECK(!memcmp(codeword, before, len),
                "rejected %zu+%zu codeword stays unchanged", p->msg_len, p->parity_len);
            continue;
        }

        /* More errors than correctable must be rejected. With two parity
         * symbols only, the syndrome pair of a double error usually matches
         * a single error somewhere else, so the decoder reports a correction
         * instead - libfec did exactly the same. What it may never report is
         * a correction that is not even a codeword of this code. */
        CHECK(p->parity_len < 3,
            "%zu errors of a %zu+%zu codeword rejected, got %d",
            errors, p->msg_len, p->parity_len, fixed);
        CHECK(is_codeword(codeword, p->msg_len, p->parity_len),
            "%d reported corrections of a %zu+%zu codeword form a codeword",
            fixed, p->msg_len, p->parity_len);
    }
}

static void test_single_parity_symbol(void)
{
    const struct rs_params p = {100, 1};
    uint8_t codeword[SCF_RS_MAX_LEN];
    uint8_t clean[SCF_RS_MAX_LEN];
    uint8_t before[SCF_RS_MAX_LEN];
    size_t len = p.msg_len + p.parity_len;

    make_codeword(codeword, &p);
    memcpy(clean, codeword, len);

    /* One parity symbol detects a single error, but cannot correct it */
    for (size_t i = 0; i < len; i++) {
        memcpy(codeword, clean, len);
        codeword[i] ^= nonzero_delta();
        memcpy(before, codeword, len);

        CHECK(scf_rs_decode(codeword, p.msg_len, p.parity_len) == -1,
            "single error at symbol %zu of a %zu+%zu codeword rejected",
            i, p.msg_len, p.parity_len);
        CHECK(!memcmp(codeword, before, len),
            "rejected %zu+%zu codeword stays unchanged", p.msg_len, p.parity_len);
    }
}

static void run(const char *name, void (*test)(void))
{
    unsigned int failures_before = failures;

    test();

    printf("%-26s %s\n", name, failures == failures_before ? "ok" : "FAILED");
}

int main(void)
{
    run("parameter checks", test_parameters);
    run("known parity", test_known_parity);
    run("encode and decode", test_roundtrip);
    run("single error per symbol", test_single_errors);
    run("correctable errors", test_correctable_errors);
    run("uncorrectable errors", test_uncorrectable_errors);
    run("single parity symbol", test_single_parity_symbol);

    printf("\n%u checks, %u failures\n", checks, failures);

    return failures != 0;
}
