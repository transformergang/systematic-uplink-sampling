/* Portable integer-only building blocks for the PAM node decision rules (no RIOT, no float). */
#ifndef PAM_COMMON_H
#define PAM_COMMON_H
#include <stdint.h>

#define PAM_NLEVEL 64
#define PAM_Q 32                      /* fixed-point fraction bits of price / token state (see README: deviation from Q16.16) */
#define PAM_ONE_Q ((int64_t)1 << PAM_Q)

typedef struct {
    const uint32_t *lut;              /* p16 LUT [nbin][64], values 1..65536 (offline, from the probability rule + round(p*65536) + clip) */
    const uint8_t  *hot;              /* per-bin flag q_j >= q_hot (EXCL online block onset) */
    uint32_t nbin;
    uint32_t p16min;                  /* round(pmin*65536): p used when token bucket < 1 */
    uint32_t G;                       /* EXCL: min non-hot units before a hot onset forces a block */
    int64_t  eu_q32;                  /* eta/step_logmu in level units, Q.32 (per upload) */
    int64_t  eb_q32;                  /* eta*b/step_logmu, Q.32 (per unit budget credit) */
    int64_t  b_q32;                   /* b, Q.32 */
    int64_t  depth_q32;               /* bucket depth, Q.32 */
    int64_t  tok0_q32;                /* initial tokens, Q.32 */
    int64_t  x0_q32;                  /* initial (logmu-lv0)/step in level units, Q.32 */
    uint32_t duty_m;                  /* DUTY_CYCLE period (units) */
    uint32_t thr_bin;                 /* THRESH: upload iff bin >= thr_bin */
} pam_cfg_t;

/* nearest level (round half up with a 2^-20-level tie bias so that mathematically exact ties, e.g. up=1,len=20,b=0.05,
 * resolve upward despite Q.32 rounding noise ~1e-9), clipped to 0..63; x in level units Q.32 */
#define PAM_TIE_EPS ((int64_t)1 << (PAM_Q - 20))
static inline uint32_t pam_level(int64_t x)
{
    int64_t k = (x + ((int64_t)1 << (PAM_Q - 1)) + PAM_TIE_EPS) >> PAM_Q;
    return (uint32_t)(k < 0 ? 0 : (k > PAM_NLEVEL - 1 ? PAM_NLEVEL - 1 : k));
}

/* Marsaglia xorshift32 13/17/5 (same recurrence as the python reference xorshift32, 32-bit masked) */
static inline uint32_t pam_xs32(uint32_t x)
{
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return x;
}

uint32_t pam_crc32(const uint8_t *buf, uint32_t len);   /* zlib CRC32 (poly 0xEDB88320) */
#endif
