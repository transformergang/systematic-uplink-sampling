#include "score.h"

void score_init(score_state_t *s)
{
    for (int i = 0; i < 2; i++) s->x1[i] = s->x2[i] = s->y1[i] = s->y2[i] = 0;
    s->base_q8 = 0; s->start = 0; s->have_base = 0;
}

static inline int32_t log2_q8(uint64_t E)
{
    if (E == 0) return 0;
    uint32_t hi = (uint32_t)(E >> 32), n;
    n = hi ? 63u - (uint32_t)__builtin_clz(hi) : 31u - (uint32_t)__builtin_clz((uint32_t)E);
    uint32_t m = n >= 6 ? (uint32_t)(E >> (n - 6)) & 63u : (uint32_t)(E << (6 - n)) & 63u;
    return (int32_t)((n << 8) + sc_loglut[m]);
}

int32_t score_unit_energy(score_state_t *s, uint32_t sh)
{
    const int32_t a0 = sc_coef[0][0], a1 = sc_coef[0][1], a2 = sc_coef[0][2], a3 = sc_coef[0][3], a4 = sc_coef[0][4];
    const int32_t c0 = sc_coef[1][0], c1 = sc_coef[1][1], c2 = sc_coef[1][2], c3 = sc_coef[1][3], c4 = sc_coef[1][4];
    int32_t x1 = s->x1[0], x2 = s->x2[0], y1 = s->y1[0], y2 = s->y2[0];
    int32_t u1 = s->x1[1], u2 = s->x2[1], v1 = s->y1[1], v2 = s->y2[1];
    uint32_t idx = s->start;
    uint64_t acc = 0;
    for (uint32_t i = 0; i < SC_UNIT_N; i++) {
        int32_t x = ((int32_t)sc_audio[idx] >> sh) << 8;   /* all filter signals in Q8 (extra 8 fractional bits suppress truncation limit cycles of the 25 Hz poles) */
        idx = (idx + 1u) & (SC_NTAB - 1u);
        int64_t t = (int64_t)a0 * x + (int64_t)a1 * x1 + (int64_t)a2 * x2 + (int64_t)a3 * y1 + (int64_t)a4 * y2;
        int32_t y = (int32_t)(t >> SC_QC);
        x2 = x1; x1 = x; y2 = y1; y1 = y;
        int64_t r = (int64_t)c0 * y + (int64_t)c1 * u1 + (int64_t)c2 * u2 + (int64_t)c3 * v1 + (int64_t)c4 * v2;
        int32_t zq = (int32_t)(r >> SC_QC);                /* Q8 */
        u2 = u1; u1 = y; v2 = v1; v1 = zq;
        int32_t z = (zq + 128) >> 8;
        acc += (uint32_t)z * (uint32_t)z;
    }
    s->x1[0] = x1; s->x2[0] = x2; s->y1[0] = y1; s->y2[0] = y2;
    s->x1[1] = u1; s->x2[1] = u2; s->y1[1] = v1; s->y2[1] = v2;
    s->start = (s->start + SC_UNIT_N) & (SC_NTAB - 1u);
    return log2_q8(acc);
}

uint32_t score_bin(score_state_t *s, int32_t L)
{
    if (!s->have_base) { s->base_q8 = L; s->have_base = 1; }
    int32_t sc = L - s->base_q8;
    s->base_q8 += (L - s->base_q8) >> 3;        /* causal: score uses the PREVIOUS baseline */
    uint32_t lo = 0, hi = 19;                    /* bin = #{thr <= sc}, binary search over 19 thresholds */
    while (lo < hi) { uint32_t mid = (lo + hi) >> 1; if (sc_thr[mid] <= sc) lo = mid + 1; else hi = mid; }
    return lo;
}

void sketch_init(sketch_t *k)
{
    uint16_t *p = &k->cnt[0][0][0];
    for (uint32_t i = 0; i < (uint32_t)SK_P * SC_NBIN * 64; i++) p[i] = 0;
    k->per = 0; k->in_per = 0;
}

void sketch_update(sketch_t *k, uint32_t bin, uint32_t level)
{
    uint16_t *p = &k->cnt[k->per][bin][level];
    if (*p != 0xFFFFu) (*p)++;
    if (++k->in_per >= SK_UPP) { k->in_per = 0; k->per = (k->per + 1u == SK_P) ? 0u : k->per + 1u; }
}
