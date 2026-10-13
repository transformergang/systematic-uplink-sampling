/* Madow plug-in decision-level benchmark (RIOT): EXCL decision vs external-p unit (RS / AM LUT) vs external-p + integer Madow step.
 * Decision only: no score, no sketch.  Structure: run each decision rule over the stored test vectors and time blocks of decisions.
 * Boards: iotlab-m3 (Cortex-M3, DWT->CYCCNT) and samr21-xpro (Cortex-M0+, SysTick 24-bit free-running down counter).
 * Every pass = VEC_N (4096) decisions over vec_bin / vec_gap; each decision is a cross-TU (non-inlined) call; the pass is timed in
 * blocks of BLK (1024) decisions so that one CYC() span stays < 2^23 cycles (M0+ 24-bit counter); blocks are summed.
 * Cycle accounting: `cyc` in PHASE END = sum over blocks of CYC_DIFF(start,end), i.e. it INCLUDES the block-function call and the
 * two CYC() reads (about `ovh` cycles per block = 4 per pass), NOT corrected.  `ovh` is printed separately (cycles of one CYC() pair);
 * the analysis may subtract (passes * 4 * ovh) if wanted.  decisions = passes * VEC_N (total decisions in the phase). */
#include <stdio.h>
#include <stdint.h>
#include "cpu.h"
#include "periph_conf.h"
#if defined(PAM_USE_ZTIMER)
#include "ztimer.h"
#define NOW_MS()      ((uint32_t)ztimer_now(ZTIMER_MSEC))
#define SLEEP_MS(ms)  ztimer_sleep(ZTIMER_MSEC, (ms))
#else
#include "xtimer.h"
#define NOW_MS()      ((uint32_t)(xtimer_now_usec64() / 1000))
#define SLEEP_MS(ms)  xtimer_usleep((uint32_t)(ms) * 1000u)
#endif
#include "vectors.h"
#include "bitexact.h"
#include "vectors_madow.h"

#define CFG          (&vec_cfg[0])          /* EXCL config: whale kappa=0.95, b=0.05 */
#define CFG_RS       (&vm_cfg[0])           /* vm_cfg index = bi*2 + marg: b=0.05 RS */
#define CFG_AM       (&vm_cfg[1])           /* b=0.05 AM */
#ifndef PHASE_S
#define PHASE_S      30
#endif
#ifndef IDLE_S
#define IDLE_S       10
#endif
#ifndef REPEATS
#define REPEATS      3
#endif
#define BLK          1024u                  /* decisions per timed block */
#define SPAN_LIMIT   (1u << 23)

#if defined(HOST_STUB)                      /* host check only: fake cycle counter from cpu.h (host_cyc), optionally masked to 24 bit */
#if defined(HOST_M0P)
#define M0P 1
static void cyc_enable(void) { }
static inline uint32_t CYC(void) { return host_cyc() & 0x00FFFFFFu; }
#define CYC_DIFF(a, b) (((b) - (a)) & 0x00FFFFFFu)
#else
#define M0P 0
static void cyc_enable(void) { }
static inline uint32_t CYC(void) { return host_cyc(); }
#define CYC_DIFF(a, b) ((uint32_t)((b) - (a)))
#endif
#elif defined(__CORTEX_M) && (__CORTEX_M == 0)
#define M0P 1
static void cyc_enable(void) { SysTick->LOAD = 0x00FFFFFFu; SysTick->VAL = 0; SysTick->CTRL = 5u; }   /* enable, core clock, no irq */
static inline uint32_t CYC(void) { return (0x00FFFFFFu - SysTick->VAL) & 0x00FFFFFFu; }
#define CYC_DIFF(a, b) (((b) - (a)) & 0x00FFFFFFu)                                                      /* valid for spans < 16.7M cycles */
#else
#define M0P 0
static void cyc_enable(void)
{
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk; DWT->CYCCNT = 0; DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
}
static inline uint32_t CYC(void) { return DWT->CYCCNT; }
#define CYC_DIFF(a, b) ((uint32_t)((b) - (a)))
#endif

static uint8_t mask[(VEC_N + 7) / 8];
static volatile uint32_t sink;
static uint32_t ovh;                    /* cost of one CYC() read pair (cycles) */

static void pu64(uint64_t v)           /* newlib-nano printf has no %llu */
{
    char b[24]; int i = 23; b[i] = 0;
    do { b[--i] = (char)('0' + (uint32_t)(v % 10u)); v /= 10u; } while (v);
    fputs(&b[i], stdout);
}
static void kv(const char *k, uint64_t v) { printf(" %s=", k); pu64(v); }

/* ---- per-pass state of all decision kinds (only the ones used by the running phase are initialised) ---- */
typedef struct { excl_state_t es; extp_state_t xs; madow_state_t ms; uint32_t busy_acc; } ctx_t;

#define NOINL __attribute__((noinline))

/* Each blk_* processes decisions [s, s+n) of the stream and returns the number of selected units (P: 0).  Cross-TU calls
 * (excl_step / extp_step / madow_step live in their own .c files) are not inlined. */
static NOINL uint32_t blk_busy(ctx_t *x, uint32_t s, uint32_t n)      
{
    uint32_t acc = x->busy_acc, nup = 0, t;
    for (t = s; t < s + n; t++) { acc += vec_bin[t]; if (acc & 0x80000000u) nup++; }
    x->busy_acc = acc; sink = acc; return nup;
}
static NOINL uint32_t blk_excl(ctx_t *x, uint32_t s, uint32_t n)         /* online forced-break path */
{
    const pam_cfg_t *c = CFG; uint32_t nup = 0, t;
    for (t = s; t < s + n; t++)
        if (excl_step(&x->es, c, vec_bin[t], excl_forced_online(&x->es, c, vec_bin[t], vec_gap[t]))) { mask[t >> 3] |= 1u << (t & 7); nup++; }
    return nup;
}
static NOINL uint32_t blk_p(ctx_t *x, const extp_cfg_t *c, uint32_t s, uint32_t n)   /* external-p unit only */
{
    uint32_t acc = 0, t;
    for (t = s; t < s + n; t++) acc += extp_step(&x->xs, c, vec_bin[t]);
    sink = acc; return 0;
}
static NOINL uint32_t blk_madow(ctx_t *x, const extp_cfg_t *c, uint32_t s, uint32_t n)  /* external-p unit + Madow step */
{
    uint32_t nup = 0, t;
    for (t = s; t < s + n; t++) {
        uint32_t p = extp_step(&x->xs, c, vec_bin[t]);
        if (madow_step(&x->ms, p, vec_gap[t])) nup++;
    }
    return nup;
}

enum { K_BUSY, K_EXCL, K_P_RS, K_P_AM, K_MADOW_RS, K_MADOW_AM };

typedef struct { uint64_t cyc, ups; uint32_t span_warn; } acc_t;

static uint32_t run_block(int kind, ctx_t *x, uint32_t s)
{
    switch (kind) {
    case K_BUSY:     return blk_busy(x, s, BLK);
    case K_EXCL:     return blk_excl(x, s, BLK);
    case K_P_RS:     return blk_p(x, CFG_RS, s, BLK);
    case K_P_AM:     return blk_p(x, CFG_AM, s, BLK);
    case K_MADOW_RS: return blk_madow(x, CFG_RS, s, BLK);
    default:         return blk_madow(x, CFG_AM, s, BLK);
    }
}

static void pass_init(int kind, ctx_t *x, uint32_t pass_n)
{
    uint32_t i;
    x->busy_acc = 0;
    switch (kind) {
    case K_EXCL:
        for (i = 0; i < sizeof mask; i++) mask[i] = 0;
        excl_init(&x->es, CFG, 0x2545F491u + pass_n); break;
    case K_P_RS: case K_P_AM: extp_init(&x->xs); break;
    case K_MADOW_RS: case K_MADOW_AM: extp_init(&x->xs); madow_init(&x->ms, 0x2545F491u + pass_n); break;
    default: break;
    }
}

static const char *kname[] = {"BUSY_BASELINE", "EXCL", "P_RS", "P_AM", "MADOW_RS", "MADOW_AM"};

static void phase(int kind, uint32_t secs)
{
    acc_t a = {0, 0, 0}; uint32_t pass_n = 0; ctx_t x;
    uint32_t t0 = NOW_MS();
    printf("PHASE %s START t=%lu\n", kname[kind], (unsigned long)t0);
    while (NOW_MS() - t0 < secs * 1000u) {
        pass_init(kind, &x, pass_n);
        for (uint32_t s = 0; s < VEC_N; s += BLK) {
            uint32_t c0 = CYC();
            uint32_t r = run_block(kind, &x, s);
            uint32_t c1 = CYC();
            uint32_t d = CYC_DIFF(c0, c1);
            if (d >= SPAN_LIMIT) { if (a.span_warn < 8) printf("CYC_SPAN_WARN phase=%s pass=%lu start=%lu span=%lu\n", kname[kind], (unsigned long)pass_n, (unsigned long)s, (unsigned long)d); a.span_warn++; }
            a.cyc += d; a.ups += r;
        }
        pass_n++;
    }
    uint32_t t1 = NOW_MS();
    printf("PHASE %s END", kname[kind]);
    kv("t", t1); kv("passes", pass_n); kv("decisions", (uint64_t)pass_n * VEC_N); kv("cyc", a.cyc); kv("ups", a.ups); kv("ovh", ovh);
    if (a.span_warn) kv("span_warn", a.span_warn);
    printf("\n");
}

static void idle(const char *label, uint32_t secs)
{
    printf("PHASE %s START t=%lu\n", label, (unsigned long)NOW_MS());
    SLEEP_MS(secs * 1000u);
    printf("PHASE %s END t=%lu\n", label, (unsigned long)NOW_MS());
}

/* ---- boot-time bit-exactness: 34 EXCL-family vectors ---- */
static int bitexact_boot(void)
{
    static const char *nm[] = {"INDEP", "EXCL", "EXCL_ONLINE", "DUTY", "THRESH"};
    int fails = 0;
    for (int ci = 0; ci < VEC_NCFG; ci++)
        for (int rule = 0; rule < 5; rule++) {
            int ns = (rule == BE_THRESH) ? 1 : VEC_NSEED;
            for (int si = 0; si < ns; si++) {
                const vec_exp_t *e = rule == BE_INDEP ? &vec_exp_indep[ci * VEC_NSEED + si] :
                    (rule == BE_DUTY ? &vec_exp_duty[ci * VEC_NSEED + si] :
                    (rule == BE_THRESH ? &vec_exp_thr[ci] : &vec_exp_excl[ci * VEC_NSEED + si]));
                be_res_t r; be_run(rule, &vec_cfg[ci], e->state, e->phase, &r);
                int ok = be_match(rule, &r, e); fails += !ok;
                printf("BITEXACT cfg=%d rule=%s seed=%d n_up=%lu crc=%08lx expect_crc=%08lx %s\n", ci, nm[rule], si,
                       (unsigned long)r.n_up, (unsigned long)r.crc, (unsigned long)e->crc, ok ? "OK" : "MISMATCH");
            }
        }
    printf("BITEXACT SUMMARY %s\n", fails ? "FAIL" : "PASS");
    return fails;
}

/* ---- boot-time bit-exactness: 32 Madow vectors, zero tolerance (n_up, mask CRC, p16 CRC, final S, final C) ---- */
static uint32_t crc_upd(uint32_t crc, uint8_t b)      /* streaming zlib CRC32 (no 16 KB p16 buffer: SAMR21 has 32 KB RAM); crc starts 0xFFFFFFFF, final ~ */
{
    crc ^= b;
    for (int k = 0; k < 8; k++) crc = (crc >> 1) ^ (0xEDB88320u & (uint32_t)(-(int32_t)(crc & 1u)));
    return crc;
}
static int bitexact_madow_boot(void)
{
    static const char *mn[] = {"RS", "AM"};
    int fails = 0;
    for (int ci = 0; ci < VM_NCFG; ci++)
        for (int si = 0; si < VM_NSEED; si++) {
            const vm_exp_t *e = &vm_exp[ci * VM_NSEED + si];
            extp_state_t xs; madow_state_t ms; uint32_t nup = 0, pcrc = 0xFFFFFFFFu, t, i;
            for (i = 0; i < sizeof mask; i++) mask[i] = 0;
            extp_init(&xs); madow_init(&ms, e->rng);
            for (t = 0; t < VEC_N; t++) {
                uint32_t p = extp_step(&xs, &vm_cfg[ci], vec_bin[t]);
                pcrc = crc_upd(pcrc, (uint8_t)p); pcrc = crc_upd(pcrc, (uint8_t)(p >> 8));
                pcrc = crc_upd(pcrc, (uint8_t)(p >> 16)); pcrc = crc_upd(pcrc, (uint8_t)(p >> 24));
                if (madow_step(&ms, p, (uint8_t)(t == 0 || vec_gap[t]))) { mask[t >> 3] |= (uint8_t)(1u << (t & 7)); nup++; }
            }
            pcrc = ~pcrc;
            uint32_t mcrc = pam_crc32(mask, sizeof mask);
            int ok = nup == e->n_up && mcrc == e->mask_crc && pcrc == e->p16_crc && xs.S == e->S && ms.C == e->C;
            fails += !ok;
            printf("BITEXACT_MADOW cfg=%d rule=%s seed=%d n_up=%lu crc=%08lx expect_crc=%08lx p16crc=%08lx expect_p16crc=%08lx", ci, mn[ci & 1], si,
                   (unsigned long)nup, (unsigned long)mcrc, (unsigned long)e->mask_crc, (unsigned long)pcrc, (unsigned long)e->p16_crc);
            printf(" S="); pu64(xs.S); printf(" expect_S="); pu64(e->S);
            printf(" C=%lu expect_C=%lu expect_n_up=%lu %s\n", (unsigned long)ms.C, (unsigned long)e->C, (unsigned long)e->n_up, ok ? "OK" : "MISMATCH");
        }
    printf("BITEXACT_MADOW SUMMARY %s (%d vectors)\n", fails ? "FAIL" : "PASS", VM_NVEC);
    return fails;
}

#ifndef SYNC_CYCLES
#define SYNC_CYCLES   5
#endif
#ifndef SYNC_BUSY_MS
#define SYNC_BUSY_MS  2000u
#endif
#ifndef SYNC_SLEEP_MS
#define SYNC_SLEEP_MS 2000u
#endif
#ifndef SYNC_TAIL_MS
#define SYNC_TAIL_MS  10000u
#endif
static void sync_pattern(void)
{
    ctx_t x; x.busy_acc = 0;
    printf("SYNC START t=%lu\n", (unsigned long)NOW_MS());
    for (int k = 0; k < SYNC_CYCLES; k++) {
        uint32_t t0 = NOW_MS();
        printf("SYNC EDGE i=%d state=busy t=%lu\n", 2 * k, (unsigned long)t0);
        while (NOW_MS() - t0 < SYNC_BUSY_MS) (void)blk_busy(&x, 0, VEC_N);
        printf("SYNC EDGE i=%d state=sleep t=%lu\n", 2 * k + 1, (unsigned long)NOW_MS());
        SLEEP_MS(SYNC_SLEEP_MS);
    }
    SLEEP_MS(SYNC_TAIL_MS);
    printf("SYNC END t=%lu\n", (unsigned long)NOW_MS());
}

int main(void)
{
    cyc_enable();
    { uint32_t s = 0; for (int i = 0; i < 1000; i++) { uint32_t a = CYC(); uint32_t b = CYC(); s += CYC_DIFF(a, b); } ovh = (s + 500) / 1000; }
    printf("MADOW_BENCH boot m0p=%d cpu_hz=%lu vec_n=%d phase_s=%d idle_s=%d repeats=%d ovh=%lu\n", M0P,
           (unsigned long)CLOCK_CORECLOCK, VEC_N, PHASE_S, IDLE_S, REPEATS, (unsigned long)ovh);
    sync_pattern();
    bitexact_boot();
    bitexact_madow_boot();
    for (int rep = 0; rep < REPEATS; rep++) {
        printf("SEQUENCE %d START t=%lu\n", rep, (unsigned long)NOW_MS());
        idle("IDLE_SLEEP", IDLE_S);
        for (int k = K_BUSY; k <= K_MADOW_AM; k++) { phase(k, PHASE_S); idle("IDLE_SLEEP", IDLE_S); }
        printf("SEQUENCE %d END t=%lu\n", rep, (unsigned long)NOW_MS());
    }
    printf("MADOW_BENCH DONE\n");
#if defined(HOST_STUB)
    return 0;                           /* host check: exit instead of idling forever */
#endif
    while (1) SLEEP_MS(1000);
    return 0;
}
