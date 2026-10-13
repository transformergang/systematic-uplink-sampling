#ifndef BASELINES_H
#define BASELINES_H
#include "pam_common.h"
typedef struct { uint32_t ctr; } duty_state_t;
void duty_init(duty_state_t *s, const pam_cfg_t *c, uint32_t phase);   /* phase in [0,m) */
int  duty_step(duty_state_t *s, const pam_cfg_t *c);                   /* upload iff (t-phase)%m==0 */
int  thresh_step(const pam_cfg_t *c, uint8_t bin);                     /* upload iff bin >= thr_bin */
#endif
