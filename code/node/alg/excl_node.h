#ifndef EXCL_NODE_H
#define EXCL_NODE_H
#include "pam_common.h"
typedef struct {
    uint32_t rng, A, u, up, blen, p16block_unused;
    int64_t x, tokens;
    uint32_t plow, level, overdraft, n_blocks;
    uint32_t nh, first;                      /* online forced-break tracker */
} excl_state_t;
void excl_init(excl_state_t *s, const pam_cfg_t *c, uint32_t rng_state);
/* forced = static block-start flag (t=0, recording gap, hot-run onset after >=G non-hot units) */
int  excl_step(excl_state_t *s, const pam_cfg_t *c, uint8_t bin, uint8_t forced);
/* compute the forced flag online from bin (hot LUT) and gap; call once per unit BEFORE excl_step */
uint8_t excl_forced_online(excl_state_t *s, const pam_cfg_t *c, uint8_t bin, uint8_t gap);
#endif
