/* External open-loop marginal probability unit (RS / AM share the code; they differ only in the weight LUT).
 *   p16_t = clip( floor( KAPQ * W[bin_t] * (K0 + t) / (C0 + S_t) ), p16min, 65536 ),  S_t = sum_{u<t} W[bin_u]  (PAST weights only)
 * then S += W[bin_t]; t += 1.  Integer only: one uint64/uint64 floor division, no float, no 128-bit.
 * Overflow contract (asserted by tools/make_vectors_madow.py): W <= 2^16, KAPQ <= 2^17, K0 + t < 2^31
 *   => KAPQ*W*(K0+t) < 2^17 * 2^16 * 2^31 = 2^64.  S <= 2^16 * t < 2^47 (uint64), C0 + S fits uint64.
 * t is a uint32 stream counter from 0; it is NOT reset at recording gaps (on-chip stream). */
#ifndef EXTP_NODE_H
#define EXTP_NODE_H
#include <stdint.h>
#define EXTP_NBIN 20
typedef struct {
    const uint32_t *w;          /* weight LUT W[bin], 1..65536 */
    uint32_t kapq;              /* round(65536 * kappa * m_ref) */
    uint32_t c0;                /* round(K0 * m_ref / sc),  sc = max_j w_j / 65536 */
    uint32_t k0;                /* prior pseudo-count */
    uint32_t p16min;            /* round(65536 * pmin) */
} extp_cfg_t;
typedef struct { uint32_t t; uint64_t S; } extp_state_t;
void extp_init(extp_state_t *s);
uint32_t extp_step(extp_state_t *s, const extp_cfg_t *c, uint8_t bin);   /* returns p16 in [p16min, 65536] */
#endif
