#ifndef SCORE_H
#define SCORE_H
/* band-level score for the vessel data, integer only: per-sample 4th-order IIR band-pass 25-500 Hz (2 Q28 biquads, DF1, int64 MAC) at SC_FS Hz,
 * sum of squares over a unit (SC_UNIT_N samples), log2 via 64-entry LUT (Q8), causal EMA baseline, bin = #thresholds <= (L - baseline).
 * The 66k whale detector is OUT OF SCOPE. Input is a stored synthetic audio table (read like an ADC DMA buffer). */
#include <stdint.h>
#include "tables.h"
typedef struct { int32_t x1[2], x2[2], y1[2], y2[2]; int32_t base_q8; uint32_t start, have_base; } score_state_t;
typedef struct { uint32_t per, in_per; uint16_t cnt[SK_P][SC_NBIN][64]; } sketch_t;
void score_init(score_state_t *s);
/* process one unit of SC_UNIT_N samples (input >> sh = per-unit gain step); returns Q8 log2 energy */
int32_t score_unit_energy(score_state_t *s, uint32_t sh);
uint32_t score_bin(score_state_t *s, int32_t L_q8);         /* baseline EMA + threshold table -> bin 0..19 */
void sketch_init(sketch_t *k);
void sketch_update(sketch_t *k, uint32_t bin, uint32_t level);
#endif
