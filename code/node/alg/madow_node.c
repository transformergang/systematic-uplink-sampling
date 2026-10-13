#include "madow_node.h"

void madow_init(madow_state_t *s, uint32_t rng_seed) { s->rng = rng_seed; s->C = 0; s->u = 0; s->first = 1; }

int madow_step(madow_state_t *s, uint32_t p16, uint8_t seg_start)
{
    if (seg_start || s->first) {
        s->first = 0;
        s->C = 0;
        s->rng = pam_xs32(s->rng);
        s->u = (s->rng >> 16) + 1u;
    }
    uint32_t Cn = s->C + p16;
    int xi = ((Cn + 65536u - s->u) >> 16) > ((s->C + 65536u - s->u) >> 16);
    s->C = Cn;
    return xi;
}
