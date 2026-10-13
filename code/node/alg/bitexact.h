/* Shared bit-exactness driver (portable; used by native test and by firmware boot pass). Needs vectors.h included first. */
#ifndef BITEXACT_H
#define BITEXACT_H
#include "indep_node.h"
#include "excl_node.h"
#include "baselines.h"
typedef struct { uint32_t n_up, crc, level, overdraft, n_blocks; int64_t tok_q32; } be_res_t;
enum { BE_INDEP, BE_EXCL, BE_EXCL_ONLINE, BE_DUTY, BE_THRESH };

static inline void be_run(int rule, const pam_cfg_t *c, uint32_t state, uint32_t phase, be_res_t *r)
{
    static uint8_t mask[(VEC_N + 7) / 8];
    for (uint32_t i = 0; i < sizeof mask; i++) mask[i] = 0;
    uint32_t nup = 0;
    indep_state_t ps; excl_state_t es; duty_state_t ds;
    indep_init(&ps, c, state); excl_init(&es, c, state); duty_init(&ds, c, phase);
    for (uint32_t t = 0; t < VEC_N; t++) {
        int xi = 0;
        switch (rule) {
        case BE_INDEP: xi = indep_step(&ps, c, vec_bin[t]); break;
        case BE_EXCL: xi = excl_step(&es, c, vec_bin[t], vec_forced[t]); break;
        case BE_EXCL_ONLINE: xi = excl_step(&es, c, vec_bin[t], excl_forced_online(&es, c, vec_bin[t], vec_gap[t])); break;
        case BE_DUTY: xi = duty_step(&ds, c); break;
        default: xi = thresh_step(c, vec_bin[t]); break;
        }
        if (xi) { mask[t >> 3] |= (uint8_t)(1u << (t & 7)); nup++; }
    }
    r->n_up = nup; r->crc = pam_crc32(mask, sizeof mask); r->level = 0; r->overdraft = 0; r->n_blocks = 0; r->tok_q32 = 0;
    if (rule == BE_INDEP) { r->level = pam_level(ps.x); r->overdraft = ps.overdraft; r->tok_q32 = ps.tokens; }
    if (rule == BE_EXCL || rule == BE_EXCL_ONLINE) {
        r->level = es.level; r->overdraft = es.overdraft; r->n_blocks = es.n_blocks; r->tok_q32 = es.tokens;
    }
}
#define BE_TOK_TOL ((int64_t)1 << 12)    /* 2^-20 tokens: float reference vs integer state */
static inline int be_match(int rule, const be_res_t *r, const vec_exp_t *e)
{
    int ok = r->n_up == e->n_up && r->crc == e->crc;
    if (rule == BE_INDEP || rule == BE_EXCL || rule == BE_EXCL_ONLINE) {
        int64_t d = r->tok_q32 - e->tok_q32; if (d < 0) d = -d;
        ok = ok && r->level == e->level && r->overdraft == e->overdraft && d <= BE_TOK_TOL;
    }
    if (rule == BE_EXCL || rule == BE_EXCL_ONLINE) ok = ok && r->n_blocks == e->n_blocks;
    return ok;
}
#endif
