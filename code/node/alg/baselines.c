#include "baselines.h"

void duty_init(duty_state_t *s, const pam_cfg_t *c, uint32_t phase) { s->ctr = (c->duty_m - phase % c->duty_m) % c->duty_m; }

int duty_step(duty_state_t *s, const pam_cfg_t *c)
{
    /* ctr = (t - phase) mod m, starts at (m-phase)%m; unit t uploads iff ctr==0 */
    int xi = (s->ctr == 0);
    s->ctr = (s->ctr + 1 >= c->duty_m) ? 0 : s->ctr + 1;
    return xi;
}

int thresh_step(const pam_cfg_t *c, uint8_t bin) { return bin >= c->thr_bin; }

uint32_t pam_crc32(const uint8_t *buf, uint32_t len)
{
    uint32_t crc = 0xFFFFFFFFu;
    for (uint32_t i = 0; i < len; i++) {
        crc ^= buf[i];
        for (int k = 0; k < 8; k++) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
    }
    return ~crc;
}
