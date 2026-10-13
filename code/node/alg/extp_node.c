#include "extp_node.h"

void extp_init(extp_state_t *s) { s->t = 0; s->S = 0; }

uint32_t extp_step(extp_state_t *s, const extp_cfg_t *c, uint8_t bin)
{
    uint64_t w = c->w[bin];
    uint64_t num = (uint64_t)c->kapq * w * ((uint64_t)c->k0 + s->t);   /* < 2^64, see extp_node.h */
    uint64_t q = num / ((uint64_t)c->c0 + s->S);                       /* uint64 floor division */
    uint32_t p;
    if (q < c->p16min) p = c->p16min;
    else if (q > 65536u) p = 65536u;
    else p = (uint32_t)q;
    s->S += w;
    s->t += 1;
    return p;
}
