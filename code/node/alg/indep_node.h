#ifndef INDEP_NODE_H
#define INDEP_NODE_H
#include "pam_common.h"
typedef struct { uint32_t rng; int64_t x, tokens; uint32_t overdraft; } indep_state_t;
void indep_init(indep_state_t *s, const pam_cfg_t *c, uint32_t rng_state);
/* one decision for a unit whose score bin is `bin`; returns 1 = upload (== independent draws with fixed-point probabilities, python reference simulate_independent(fixed=True)) */
int  indep_step(indep_state_t *s, const pam_cfg_t *c, uint8_t bin);
#endif
