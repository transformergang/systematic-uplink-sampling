#include "indep_node.h"

void indep_init(indep_state_t *s, const pam_cfg_t *c, uint32_t rng_state)
{
    s->rng = rng_state; s->x = c->x0_q32; s->tokens = c->tok0_q32; s->overdraft = 0;
}

int indep_step(indep_state_t *s, const pam_cfg_t *c, uint8_t bin)
{
    uint32_t p16;
    if (s->tokens >= PAM_ONE_Q) p16 = c->lut[(uint32_t)bin * PAM_NLEVEL + pam_level(s->x)];
    else                         p16 = c->p16min;
    s->rng = pam_xs32(s->rng);
    int xi = (s->rng >> 16) < p16;                       /* top16 < p16 */
    int64_t t = s->tokens + c->b_q32 - (xi ? PAM_ONE_Q : 0);
    s->tokens = t > c->depth_q32 ? c->depth_q32 : t;
    if (s->tokens < -PAM_ONE_Q) s->overdraft = 1;
    s->x += (xi ? c->eu_q32 : 0) - c->eb_q32;            /* eta*(xi-b) */
    return xi;
}
