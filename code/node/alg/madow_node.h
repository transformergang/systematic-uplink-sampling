/* Integer single-start Madow (systematic) sampling, one random start per recording segment, open loop.
 * seg_start (first unit or recording gap): C := 0; rng := xorshift32(rng); u := (rng >> 16) + 1  (1..65536).
 * every unit (p16 in 1..65536): Cn := C + p16; select iff ((Cn + 65536 - u) >> 16) > ((C + 65536 - u) >> 16); C := Cn.
 * p16 <= 65536 => at most one selection per unit.  C is uint32: a segment must have sum p16 < 2^32 - 2^16 (<= 65535 units at p=1). */
#ifndef MADOW_NODE_H
#define MADOW_NODE_H
#include "pam_common.h"
typedef struct { uint32_t rng, C, u, first; } madow_state_t;
void madow_init(madow_state_t *s, uint32_t rng_seed);
int madow_step(madow_state_t *s, uint32_t p16, uint8_t seg_start);
#endif
