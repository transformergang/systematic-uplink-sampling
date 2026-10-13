# Cost on the node

Measured on FIT IoT-LAB boards with RIOT: STM32F103 Cortex-M3 at 72 MHz (two boards; INA226 power monitor) and a SAMR21 Cortex-M0+ at 48 MHz (no power monitor, so no energy). **What was measured:** the decision rules (one call per unit, over 4096 stored score bins, three repetitions of 30 s) and the integer band-level scorer on stored audio. **What was not measured:** acquisition, storage, radio, or any detector; the whale scores of the paper were computed offline. Energy is board input power above the sleep level, so it excludes the sleep floor. All boards reproduced the fixed-point reference bit-exactly on all 66 test vectors (34 for the closed-loop family, 32 for the single-start rule). Cycle counts were identical over repetitions and over the two Cortex-M3 boards.

## Decisions

Cycles per decision include the call overhead of a non-inlined function. The rows marked "earlier" come from an earlier run on one Cortex-M3 board and were not repeated on the Cortex-M0+; the other rows ("main") come from the run on two Cortex-M3 boards and one Cortex-M0+ board. In the earlier run EXCL took 208.0 cycles.

| Rule | Cortex-M3 cycles | Cortex-M0+ cycles | Energy per decision on the M3 (uJ) | Run |
|---|---:|---:|---:|---:|
| Trigger (top scores) | 44.5 | not measured | 0.034 | earlier |
| Duty cycle | 39.4 | not measured | 0.033 | earlier |
| Closed-loop rule, independent draws (per-unit price update) | 194.1 | not measured | 0.165 | earlier |
| EXCL (closed-loop rule) | 204.5 | 232.2 | 0.185 | main |
| Independent, open-loop probabilities (Robust Sampling) | 257.4 | 762.7 | 0.238 | main |
| Independent, open-loop probabilities (Active Measurement) | 257.1 | 755.4 | 0.238 | main |
| Madow, open-loop (Robust Sampling) | 306.6 | 817.2 | 0.280 | main |
| Madow, open-loop (Active Measurement) | 306.3 | 809.9 | 0.280 | main |

## Band-level score, per 10-s unit (2 kHz, 20000 samples)

| Quantity | Value |
|---|---:|
| Cycles per unit, Cortex-M3 | 3268256 (45.4 ms) |
| Cycles per unit, Cortex-M0+ | 14154015 (294.9 ms) |
| Energy per unit on the M3, above sleep (two boards, uJ) | 3400 and 3354 |
| Share of the unit cycles taken by the EXCL decision, M3 | 0.0063 % |

## Object-code size

Cross-compiled with `-Os -std=c99` for the Cortex-M3 and Cortex-M0+ (bytes of executable code; the decision units do not use data or bss).

| Unit | Cortex-M3 (.text) | Cortex-M0+ (.text) |
|---|---:|---:|
| Open-loop probability unit (look-up, running sum, one 64-bit division) | 134 | 116 |
| Single-start step (Madow) | 96 | 84 |

Sizes of the EXCL decision unit (468 bytes on the Cortex-M3, 446 on the Cortex-M0+) are those reported with the same compilation setup.

## Uplink share of the sketch

The sketch holds the counts of all units per (period, bin, probability level) cell that is occupied. Bytes per unit assume 4 bytes per occupied cell and are measured on a replication of the deployed (fixed-point) EXCL on the test deployments. The share of the uplink is the sketch divided by the sketch plus the uploaded segments, when every uploaded segment fills one 340-byte Iridium short-burst message plus 8 bytes of metadata (time stamp, bin, probability level, period). Units per day: 1440 (whale, 60 s), 240 (vessel, 10-s clips as in the dataset), 144 (Gray's Reef, 10 min).

| Dataset | b | Occupied cells per period (mean) | Sketch bytes per unit | Sketch bytes per day | Segment bytes per day (340 + 8 B each) | Share of the uplink |
|---|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | 31.7 | 0.63 | 905 | 25056 | 3.5 % |
| Whale | 0.10 | 37.2 | 0.74 | 1067 | 50112 | 2.1 % |
| Vessel | 0.05 | 37.0 | 0.61 | 147 | 4176 | 3.4 % |
| Vessel | 0.10 | 38.8 | 0.65 | 157 | 8352 | 1.8 % |
| Gray's Reef | 0.05 | 28.7 | 0.76 | 110 | 2506 | 4.2 % |
| Gray's Reef | 0.10 | 32.2 | 0.87 | 126 | 5011 | 2.4 % |
