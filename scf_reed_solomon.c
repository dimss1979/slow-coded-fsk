/* Agent: OpenCode, model: MiMo-V2.6-Flash Free */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "scf_private.h"

/* Reed-Solomon codec for a single code:
 *  - GF(256) with the primitive polynomial 0x11d
 *  - generator polynomial roots alpha^1, alpha^2, ...
 *  - hard decision decoding only, erasures are not supported
 *
 * The message and parity lengths are flexible: any combination with at
 * least one message and one parity symbol that fits into a codeword of
 * SCF_RS_MAX_LEN symbols can be used. A codeword of msg_len + parity_len
 * symbols corrects up to parity_len / 2 symbol errors.
 */

#define RS_FIELD_SIZE 256
#define RS_MODULUS (RS_FIELD_SIZE - 1)
#define RS_POLY 0x11d
#define RS_ROOT_EXP 1 /* Exponent of the first generator root */

static uint8_t rs_exp[2 * SCF_RS_MAX_LEN];
static uint8_t rs_log[RS_FIELD_SIZE];
static bool rs_tables_ready;

static void rs_tables_init(void)
{
    if (rs_tables_ready)
        return;

    unsigned int x = 1;
    for (unsigned int i = 0; i < RS_MODULUS; i++) {
        rs_exp[i] = (uint8_t) x;
        rs_log[x] = (uint8_t) i;

        x <<= 1;
        if (x & RS_FIELD_SIZE)
            x ^= RS_POLY;
    }
    /* Duplicate the table so that no exponent needs to be wrapped */
    for (unsigned int i = RS_MODULUS; i < sizeof(rs_exp); i++)
        rs_exp[i] = rs_exp[i - RS_MODULUS];

    rs_log[0] = 0; /* Unused, keeps the table defined for zero */
    rs_tables_ready = true;
}

static uint8_t rs_mul(uint8_t a, uint8_t b)
{
    if (!a || !b)
        return 0;

    return rs_exp[rs_log[a] + rs_log[b]];
}

static uint8_t rs_div(uint8_t a, uint8_t b)
{
    if (!a)
        return 0;

    return rs_exp[rs_log[a] + RS_MODULUS - rs_log[b]];
}

/* alpha^(-e) */
static uint8_t rs_exp_neg(unsigned int e)
{
    return rs_exp[RS_MODULUS - e % RS_MODULUS];
}

/* Generator polynomial with the roots alpha^1 .. alpha^parity_len,
 * coefficients indexed by the power of x, monic of degree parity_len */
static void rs_generator(uint8_t *gen, size_t parity_len)
{
    memset(gen, 0, (parity_len + 1) * sizeof(gen[0]));
    gen[0] = 1;

    for (size_t i = 0; i < parity_len; i++) {
        uint8_t root = rs_exp[RS_ROOT_EXP + i];
        uint8_t top = gen[i];

        for (size_t j = i; j > 0; j--)
            gen[j] = gen[j - 1] ^ rs_mul(root, gen[j]);

        gen[0] = rs_mul(root, gen[0]);
        gen[i + 1] = top;
    }
}

int scf_rs_encode(uint8_t *codeword, size_t msg_len, size_t parity_len)
{
    if (!codeword || !msg_len || !parity_len || msg_len + parity_len > SCF_RS_MAX_LEN)
        return -1;

    rs_tables_init();

    uint8_t gen[SCF_RS_MAX_LEN];
    uint8_t rem[SCF_RS_MAX_LEN] = {0};

    rs_generator(gen, parity_len);

    /* Polynomial division of the message shifted by parity_len symbols,
     * the remainder becomes the parity appended to the message */
    for (size_t i = 0; i < msg_len; i++) {
        uint8_t feedback = codeword[i] ^ rem[0];

        memmove(&rem[0], &rem[1], parity_len - 1);
        rem[parity_len - 1] = 0;

        if (feedback) {
            for (size_t j = 0; j < parity_len; j++)
                rem[j] ^= rs_mul(feedback, gen[parity_len - 1 - j]);
        }
    }

    memcpy(&codeword[msg_len], rem, parity_len);

    return 0;
}

int scf_rs_decode(uint8_t *codeword, size_t msg_len, size_t parity_len)
{
    if (!codeword || !msg_len || parity_len < 1 || msg_len + parity_len > SCF_RS_MAX_LEN)
        return -1;

    rs_tables_init();

    size_t len = msg_len + parity_len;
    uint8_t synd[SCF_RS_MAX_LEN];
    bool clean = true;

    for (size_t i = 0; i < parity_len; i++) {
        uint8_t root = rs_exp[RS_ROOT_EXP + i];
        uint8_t sum = 0;

        for (size_t j = 0; j < len; j++)
            sum = rs_mul(sum, root) ^ codeword[j];

        synd[i] = sum;
        if (sum)
            clean = false;
    }

    if (clean)
        return 0;

    /* A single parity symbol detects, but cannot correct any error */
    if (parity_len < 2)
        return -1;

    /* Berlekamp-Massey: derive the error locator lambda(x) */
    uint8_t lambda[SCF_RS_MAX_LEN + 1] = {1};
    uint8_t previous[SCF_RS_MAX_LEN + 1] = {1};
    uint8_t backup[SCF_RS_MAX_LEN + 1];
    size_t errors = 0;
    size_t delay = 1;
    uint8_t scale = 1;

    for (size_t i = 0; i < parity_len; i++) {
        uint8_t disc = synd[i];

        for (size_t j = 1; j <= errors; j++)
            disc ^= rs_mul(lambda[j], synd[i - j]);

        if (!disc) {
            delay++;
            continue;
        }

        memcpy(backup, lambda, sizeof(lambda));

        uint8_t factor = rs_div(disc, scale);
        for (size_t j = 0; j + delay < sizeof(lambda); j++) {
            if (previous[j])
                lambda[j + delay] ^= rs_mul(factor, previous[j]);
        }

        if (2 * errors <= i) {
            errors = i + 1 - errors;
            memcpy(previous, backup, sizeof(previous));
            scale = disc;
            delay = 1;
        } else {
            delay++;
        }
    }

    /* More errors than this code can correct */
    if (!errors || 2 * errors > parity_len)
        return -1;

    /* Chien search: evaluate lambda at the inverse of every symbol root */
    size_t positions[SCF_RS_MAX_LEN];
    uint8_t values[SCF_RS_MAX_LEN];
    size_t found = 0;

    for (size_t p = 0; p < len; p++) {
        unsigned int symbol_root = (unsigned int) (len - 1 - p);
        uint8_t sum = 0;

        for (size_t j = 0; j <= errors; j++) {
            if (lambda[j])
                sum ^= rs_mul(lambda[j], rs_exp_neg(symbol_root * j));
        }

        if (!sum) {
            positions[found] = p;
            found++;
        }
    }

    /* Not as many error locations as the locator claims */
    if (found != errors)
        return -1;

    /* Error evaluator omega(x) = syndrome(x) * lambda(x) mod x^parity_len */
    uint8_t omega[SCF_RS_MAX_LEN] = {0};

    for (size_t i = 0; i < parity_len; i++) {
        for (size_t j = 0; j <= i; j++)
            omega[i] ^= rs_mul(synd[i - j], lambda[j]);
    }

    /* Forney algorithm: the magnitude of every located error */
    for (size_t e = 0; e < found; e++) {
        unsigned int symbol_root = (unsigned int) (len - 1 - positions[e]);
        uint8_t num = 0;
        uint8_t den = 0;

        for (size_t i = 0; i < parity_len; i++) {
            if (omega[i])
                num ^= rs_mul(omega[i], rs_exp_neg(symbol_root * i));
        }

        /* Derivative of lambda, only odd powers remain in characteristic 2 */
        for (size_t j = 1; j <= errors; j += 2) {
            if (lambda[j])
                den ^= rs_mul(lambda[j], rs_exp_neg(symbol_root * (j - 1)));
        }

        if (!den)
            return -1;

        values[e] = rs_div(num, den);
    }

    for (size_t e = 0; e < found; e++)
        codeword[positions[e]] ^= values[e];

    return (int) found;
}
