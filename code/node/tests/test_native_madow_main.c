/* Native bit-exact driver: 34 EXCL-family vectors (../vectors.h, alg/bitexact.h) + 32 Madow vectors (vectors_madow.h) + 4 self-checks. Zero tolerance on the Madow integers. */
#include <stdio.h>
#include <string.h>
#include "vectors.h"
#include "bitexact.h"
#include "vectors_madow.h"

static uint32_t crc_p16(const uint32_t *p, uint32_t n)
{
    static uint8_t buf[4 * VEC_N];
    for (uint32_t i = 0; i < n; i++) { buf[4 * i] = (uint8_t)p[i]; buf[4 * i + 1] = (uint8_t)(p[i] >> 8); buf[4 * i + 2] = (uint8_t)(p[i] >> 16); buf[4 * i + 3] = (uint8_t)(p[i] >> 24); }
    return pam_crc32(buf, 4 * n);
}

static int selfcheck_hand(void)    /* (i) */
{
    extp_state_t s; extp_init(&s); int ok = 1;
    for (int t = 0; t < VH_N; t++) { uint32_t p = extp_step(&s, &vh_cfg, vh_bin[t]); if (p != vh_p16[t]) { ok = 0; printf("  hand t=%d p16=%u exp=%u\n", t, p, vh_p16[t]); } }
    if (s.t != VH_N || s.S != VH_S) { ok = 0; printf("  hand final t=%u S=%llu\n", s.t, (unsigned long long)s.S); }
    return ok;
}

static int selfcheck_p16crc(void)  /* (ii) all 4 cfgs: p16 sequence CRC + final S == Python (also covered per vector); p16 in range */
{
    static uint32_t p[VEC_N]; int ok = 1;
    for (int ci = 0; ci < VM_NCFG; ci++) {
        extp_state_t s; extp_init(&s);
        for (int t = 0; t < VEC_N; t++) { p[t] = extp_step(&s, &vm_cfg[ci], vec_bin[t]); if (p[t] < vm_cfg[ci].p16min || p[t] > 65536u) ok = 0; }
        const vm_exp_t *e = &vm_exp[ci * VM_NSEED];
        if (crc_p16(p, VEC_N) != e->p16_crc || s.S != e->S) ok = 0;
    }
    return ok;
}

static uint32_t tr_rng;
static uint32_t tr_next(void) { tr_rng = pam_xs32(tr_rng); return tr_rng; }

static int selfcheck_madow_exact(void)   /* (iii) */
{
    int ok = 1; madow_state_t m;
    madow_init(&m, 12345u);
    for (int t = 0; t < 1000; t++) if (!madow_step(&m, 65536u, t == 0)) ok = 0;         /* p=1 everywhere: every unit selected */
    tr_rng = 987654321u;
    for (int k = 0; k < 200; k++) {
        uint32_t L = 1 + tr_next() % 300, seed = 1 + tr_next() % 0xFFFFFFFEu; uint32_t sum = 0, M, p; int cnt = 0; static uint32_t pv[400];
        for (uint32_t i = 0; i + 1 < L; i++) { uint32_t mode = tr_next() % 4; pv[i] = mode == 0 ? 65536u : (mode == 1 ? 1u : 1 + tr_next() % 65536u); sum += pv[i]; }
        p = 65536u - (sum % 65536u); if (p == 0) p = 65536u;                               /* last unit completes the integer mass */
        pv[L - 1] = p; sum += p; M = sum / 65536u;
        madow_init(&m, seed);
        for (uint32_t i = 0; i < L; i++) cnt += madow_step(&m, pv[i], i == 0);
        if ((uint32_t)cnt != M) { ok = 0; printf("  segment %d: L=%u mass M=%u selected=%d\n", k, L, M, cnt); }
    }
    return ok;
}

static int selfcheck_gap(void)           /* (iv) */
{
    madow_state_t m; int ok = 1; uint32_t u_before = 0, rng_before = 0, C_before = 0;
    madow_init(&m, 424242u);
    for (int t = 0; t < 100; t++) {
        uint8_t gap = (t == 0 || t == 50);
        if (t == 50) { u_before = m.u; rng_before = m.rng; C_before = m.C; if (C_before != 50u * 1000u) ok = 0; }
        uint32_t rng_pre = m.rng; (void)madow_step(&m, 1000u, gap);
        if (t == 50) {
            if (m.C != 1000u) ok = 0;                                              /* C cleared at the gap, then += p16 */
            if (m.rng != pam_xs32(rng_before) || m.u != (m.rng >> 16) + 1u) ok = 0; /* exactly one new draw */
            if (m.u == u_before) { printf("  note: u unchanged by redraw (coincidence)\n"); }
        } else if (m.rng != (gap ? pam_xs32(rng_pre) : rng_pre)) ok = 0;            /* no draw on non-start units */
    }
    return ok;
}

int main(void)
{
    int fails = 0, total = 0;
    static const char *nm[] = {"INDEP", "EXCL", "EXCL_ONLINE", "DUTY", "THRESH"};
    for (int ci = 0; ci < VEC_NCFG; ci++) {
        for (int rule = 0; rule < 5; rule++) {
            int ns = (rule == BE_THRESH) ? 1 : VEC_NSEED;
            for (int si = 0; si < ns; si++) {
                const vec_exp_t *e = rule == BE_INDEP ? &vec_exp_indep[ci * VEC_NSEED + si] :
                    (rule == BE_DUTY ? &vec_exp_duty[ci * VEC_NSEED + si] :
                    (rule == BE_THRESH ? &vec_exp_thr[ci] : &vec_exp_excl[ci * VEC_NSEED + si]));
                be_res_t r; be_run(rule, &vec_cfg[ci], e->state, e->phase, &r);
                int ok = be_match(rule, &r, e); total++; fails += !ok;
                printf("cfg%d %-11s seed%d n_up=%u/%u crc=%08x/%08x %s\n", ci, nm[rule], si, r.n_up, e->n_up, r.crc, e->crc, ok ? "OK" : "MISMATCH");
            }
        }
    }
    static const char *mn[] = {"RS", "AM"};
    for (int ci = 0; ci < VM_NCFG; ci++) {
        for (int si = 0; si < VM_NSEED; si++) {
            const vm_exp_t *e = &vm_exp[ci * VM_NSEED + si];
            static uint32_t p[VEC_N]; static uint8_t mask[(VEC_N + 7) / 8];
            memset(mask, 0, sizeof mask);
            extp_state_t xs; madow_state_t ms; uint32_t nup = 0;
            extp_init(&xs); madow_init(&ms, e->rng);
            for (uint32_t t = 0; t < VEC_N; t++) {
                p[t] = extp_step(&xs, &vm_cfg[ci], vec_bin[t]);
                if (madow_step(&ms, p[t], (uint8_t)(t == 0 || vec_gap[t]))) { mask[t >> 3] |= (uint8_t)(1u << (t & 7)); nup++; }
            }
            uint32_t mc = pam_crc32(mask, sizeof mask), pc = crc_p16(p, VEC_N);
            int ok = nup == e->n_up && mc == e->mask_crc && pc == e->p16_crc && xs.S == e->S && ms.C == e->C;
            total++; fails += !ok;
            printf("cfg%d MADOW_%s b=%s seed%d n_up=%u/%u maskcrc=%08x/%08x p16crc=%08x/%08x S=%llu/%llu C=%u/%u %s\n", ci, mn[ci & 1], ci < 2 ? "0.05" : "0.10", si, nup, e->n_up,
                   mc, e->mask_crc, pc, e->p16_crc, (unsigned long long)xs.S, (unsigned long long)e->S, ms.C, e->C, ok ? "OK" : "MISMATCH");
        }
    }
    int c1 = selfcheck_hand(), c2 = selfcheck_p16crc(), c3 = selfcheck_madow_exact(), c4 = selfcheck_gap();
    printf("SELFCHECK (i) extp hand example vs Python: %s\n", c1 ? "OK" : "FAIL");
    printf("SELFCHECK (ii) 4096-unit p16 CRC + final S == Python, 4 cfgs: %s\n", c2 ? "OK" : "FAIL");
    printf("SELFCHECK (iii) Madow p=1 all selected; integer-mass segments select exactly M (200 random segments): %s\n", c3 ? "OK" : "FAIL");
    printf("SELFCHECK (iv) gap in mid-stream clears C, exactly one new u: %s\n", c4 ? "OK" : "FAIL");
    int sc = c1 + c2 + c3 + c4; fails += 0;
    if (fails || sc != 4) { printf("FAIL: %d/%d vectors bit-exact (34 EXCL-family + 32 Madow), selfchecks %d/4\n", total - fails, total, sc); return 1; }
    printf("PASS: %d/%d vectors bit-exact (34 EXCL-family + 32 Madow), selfchecks 4/4\n", total - fails, total);
    return 0;
}
