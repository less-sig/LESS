/// Reference vector generator for the rust port.
///
/// Drives the C implementation from the same deterministic SHAKE based
/// `platform_csprng_state` the rust port uses, so both sides consume the very
/// same random stream and the outputs can be compared byte for byte.
///
/// Build with `kat/run.sh`, which compiles the sources of
/// `Reference_Implementation` without modifying them.
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "LESS.h"
#include "rng.h"

static void put_hex(const char *name, const unsigned char *p, size_t len) {
    printf("%s = ", name);
    for (size_t i = 0; i < len; i++) {
        printf("%02x", p[i]);
    }
    printf("\n");
}

int main(void) {
    /* the seed both implementations start from */
    unsigned char seed[48];
    for (unsigned i = 0; i < sizeof(seed); i++) {
        seed[i] = (unsigned char) i;
    }
    init_randombytes(seed, sizeof(seed));

    printf("N = %u\n", (unsigned) N);
    printf("K = %u\n", (unsigned) K);
    printf("T = %u\n", (unsigned) T);
    printf("W = %u\n", (unsigned) W);
    printf("SEED_LENGTH_BYTES = %u\n", (unsigned) SEED_LENGTH_BYTES);
    printf("RREF_MAT_PACKEDBYTES = %u\n", (unsigned) RREF_MAT_PACKEDBYTES);

    prikey_t sk;
    pubkey_t pk;
    memset(&sk, 0, sizeof(sk));
    memset(&pk, 0, sizeof(pk));

    LESS_keygen(&sk, &pk);
    put_hex("sk", (const unsigned char *) &sk, sizeof(sk));
    put_hex("pk", (const unsigned char *) &pk, sizeof(pk));

    /* the message both implementations sign */
    static const char msg[] = "LESS rust port reference vector";
    const size_t mlen = sizeof(msg) - 1;
    put_hex("msg", (const unsigned char *) msg, mlen);

    sign_t sig;
    memset(&sig, 0, sizeof(sig));
    const size_t leaves = LESS_sign(&sk, msg, mlen, &sig);

    printf("leaves = %u\n", (unsigned) leaves);
    put_hex("digest", sig.digest, HASH_DIGEST_LENGTH);
    put_hex("salt", sig.salt, HASH_DIGEST_LENGTH);
    put_hex("cf_monom_actions", (const unsigned char *) sig.cf_monom_actions, (size_t) W * N8);
    put_hex("seed_storage", sig.seed_storage, SEED_TREE_MAX_PUBLISHED_BYTES);

    const int ok = LESS_verify(&pk, msg, mlen, &sig);
    printf("verify = %d\n", ok);

    return 0;
}
