#include "excl_node.h"

void excl_init(excl_state_t *s, const pam_cfg_t *c, uint32_t rng_state)
{
    s->rng = rng_state; s->A = 0; s->u = 0; s->up = 0; s->blen = 0; s->p16block_unused = 0;
    s->x = c->x0_q32; s->tokens = c->tok0_q32;
    s->level = pam_level(s->x); s->plow = 0; s->overdraft = 0; s->n_blocks = 0;
    s->nh = 0; s->first = 1;
}

uint8_t excl_forced_online(excl_state_t *s, const pam_cfg_t *c, uint8_t bin, uint8_t gap)
{
    uint8_t hot = c->hot[bin];
    if (gap) s->nh = 0;
    uint8_t f = s->first || gap || (hot && s->nh >= c->G);
    s->first = 0;
    if (hot) s->nh = 0; else if (s->nh < 255) s->nh++;
    return f;
}

int excl_step(excl_state_t *s, const pam_cfg_t *c, uint8_t bin, uint8_t forced)
{
    const uint32_t *row = c->lut + (uint32_t)bin * PAM_NLEVEL;
    uint32_t p = s->plow ? c->p16min : row[s->level];       /* current block's frozen state */
    if (forced || s->A + p > 65536u) {                      /* open a new block: close the previous one */
        int64_t t = s->tokens + (int64_t)s->blen * c->b_q32 - ((int64_t)s->up << PAM_Q);
        s->tokens = t > c->depth_q32 ? c->depth_q32 : t;
        if (s->tokens < -PAM_ONE_Q) s->overdraft = 1;
        s->x += (int64_t)s->up * c->eu_q32 - (int64_t)s->blen * c->eb_q32;
        s->level = pam_level(s->x);
        s->plow = s->tokens < PAM_ONE_Q;
        p = s->plow ? c->p16min : row[s->level];
        s->A = 0; s->up = 0; s->blen = 0;
        s->rng = pam_xs32(s->rng);
        s->u = (s->rng >> 16) + 1u;                         /* uniform on 1..65536 */
        s->n_blocks++;
    }
    uint32_t An = s->A + p;
    int xi = (s->A < s->u) && (s->u <= An);
    s->up += (uint32_t)xi; s->blen++; s->A = An;
    return xi;
}
