# Robustness checks

## 1. Event definition

An event is a maximal run of positive units in which gaps of at most g negative units are merged (g = 1 in the paper; g = 0 and g = 2 here). Closed-loop rule and the open-loop rules (RS = Robust Sampling, AM = Active Measurement) with independent (Ind), EXCL and Madow draws; 2000 replications on the test data; event recall pooled over events.

| Dataset | g (merged gap) | Events | Median length | Mean length |
|---|---:|---:|---:|---:|
| Whale | 0 | 7022 | 2.0 | 3.18 |
| Whale | 1 | 5620 | 4.0 | 4.23 |
| Whale | 2 | 5096 | 5.0 | 4.87 |
| Vessel | 0 | 1262 | 3.0 | 3.63 |
| Vessel | 1 | 1109 | 3.0 | 4.27 |
| Vessel | 2 | 1001 | 3.0 | 4.95 |

Event recall:

| Dataset | b | g | Closed-loop Ind | Closed-loop EXCL | RS Ind | RS EXCL | RS Madow | AM Ind | AM EXCL | AM Madow |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | 0 | 0.1788 | 0.1984 | 0.1547 | 0.1693 | 0.1711 | 0.1936 | 0.2189 | 0.2210 |
| Whale | 0.05 | 1 | 0.2112 | 0.2358 | 0.1846 | 0.2037 | 0.2061 | 0.2282 | 0.2597 | 0.2621 |
| Whale | 0.05 | 2 | 0.2256 | 0.2510 | 0.1979 | 0.2177 | 0.2207 | 0.2432 | 0.2761 | 0.2790 |
| Whale | 0.1 | 0 | 0.3129 | 0.3688 | 0.2800 | 0.3243 | 0.3288 | 0.3361 | 0.4029 | 0.4112 |
| Whale | 0.1 | 1 | 0.3622 | 0.4317 | 0.3267 | 0.3818 | 0.3865 | 0.3878 | 0.4717 | 0.4807 |
| Whale | 0.1 | 2 | 0.3838 | 0.4583 | 0.3472 | 0.4059 | 0.4108 | 0.4105 | 0.5012 | 0.5110 |
| Vessel | 0.05 | 0 | 0.2043 | 0.2273 | 0.1779 | 0.1947 | 0.1998 | 0.2593 | 0.2974 | 0.3059 |
| Vessel | 0.05 | 1 | 0.2257 | 0.2538 | 0.1966 | 0.2173 | 0.2244 | 0.2852 | 0.3300 | 0.3415 |
| Vessel | 0.05 | 2 | 0.2428 | 0.2744 | 0.2122 | 0.2354 | 0.2451 | 0.3046 | 0.3525 | 0.3656 |
| Vessel | 0.1 | 0 | 0.3615 | 0.4185 | 0.3162 | 0.3611 | 0.3746 | 0.4373 | 0.4999 | 0.5308 |
| Vessel | 0.1 | 1 | 0.3913 | 0.4544 | 0.3435 | 0.3944 | 0.4113 | 0.4698 | 0.5353 | 0.5657 |
| Vessel | 0.1 | 2 | 0.4136 | 0.4776 | 0.3642 | 0.4174 | 0.4353 | 0.4924 | 0.5572 | 0.5872 |

Gain over independent draws of the same rule (difference of event recall, its standard error, and the ratio of recalls):

| Dataset | b | g | Rule | Draw | Recall difference | SE | Recall ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | 0 | closed-loop | EXCL | +0.0196 | 0.0001 | 1.110 |
| Whale | 0.05 | 0 | RS | EXCL | +0.0146 | 0.0001 | 1.095 |
| Whale | 0.05 | 0 | RS | Madow | +0.0164 | 0.0001 | 1.106 |
| Whale | 0.05 | 0 | AM | EXCL | +0.0252 | 0.0001 | 1.130 |
| Whale | 0.05 | 0 | AM | Madow | +0.0273 | 0.0001 | 1.141 |
| Whale | 0.05 | 1 | closed-loop | EXCL | +0.0245 | 0.0001 | 1.116 |
| Whale | 0.05 | 1 | RS | EXCL | +0.0191 | 0.0002 | 1.104 |
| Whale | 0.05 | 1 | RS | Madow | +0.0216 | 0.0002 | 1.117 |
| Whale | 0.05 | 1 | AM | EXCL | +0.0314 | 0.0002 | 1.138 |
| Whale | 0.05 | 1 | AM | Madow | +0.0339 | 0.0002 | 1.149 |
| Whale | 0.05 | 2 | closed-loop | EXCL | +0.0254 | 0.0001 | 1.113 |
| Whale | 0.05 | 2 | RS | EXCL | +0.0199 | 0.0002 | 1.101 |
| Whale | 0.05 | 2 | RS | Madow | +0.0229 | 0.0002 | 1.115 |
| Whale | 0.05 | 2 | AM | EXCL | +0.0329 | 0.0002 | 1.135 |
| Whale | 0.05 | 2 | AM | Madow | +0.0358 | 0.0002 | 1.147 |
| Whale | 0.1 | 0 | closed-loop | EXCL | +0.0559 | 0.0001 | 1.179 |
| Whale | 0.1 | 0 | RS | EXCL | +0.0442 | 0.0002 | 1.158 |
| Whale | 0.1 | 0 | RS | Madow | +0.0487 | 0.0002 | 1.174 |
| Whale | 0.1 | 0 | AM | EXCL | +0.0668 | 0.0001 | 1.199 |
| Whale | 0.1 | 0 | AM | Madow | +0.0751 | 0.0001 | 1.223 |
| Whale | 0.1 | 1 | closed-loop | EXCL | +0.0695 | 0.0001 | 1.192 |
| Whale | 0.1 | 1 | RS | EXCL | +0.0551 | 0.0002 | 1.169 |
| Whale | 0.1 | 1 | RS | Madow | +0.0598 | 0.0002 | 1.183 |
| Whale | 0.1 | 1 | AM | EXCL | +0.0839 | 0.0002 | 1.216 |
| Whale | 0.1 | 1 | AM | Madow | +0.0929 | 0.0002 | 1.240 |
| Whale | 0.1 | 2 | closed-loop | EXCL | +0.0745 | 0.0001 | 1.194 |
| Whale | 0.1 | 2 | RS | EXCL | +0.0587 | 0.0002 | 1.169 |
| Whale | 0.1 | 2 | RS | Madow | +0.0636 | 0.0002 | 1.183 |
| Whale | 0.1 | 2 | AM | EXCL | +0.0907 | 0.0002 | 1.221 |
| Whale | 0.1 | 2 | AM | Madow | +0.1005 | 0.0002 | 1.245 |
| Vessel | 0.05 | 0 | closed-loop | EXCL | +0.0230 | 0.0003 | 1.113 |
| Vessel | 0.05 | 0 | RS | EXCL | +0.0168 | 0.0003 | 1.094 |
| Vessel | 0.05 | 0 | RS | Madow | +0.0219 | 0.0003 | 1.123 |
| Vessel | 0.05 | 0 | AM | EXCL | +0.0381 | 0.0003 | 1.147 |
| Vessel | 0.05 | 0 | AM | Madow | +0.0466 | 0.0003 | 1.180 |
| Vessel | 0.05 | 1 | closed-loop | EXCL | +0.0281 | 0.0003 | 1.125 |
| Vessel | 0.05 | 1 | RS | EXCL | +0.0206 | 0.0004 | 1.105 |
| Vessel | 0.05 | 1 | RS | Madow | +0.0278 | 0.0003 | 1.141 |
| Vessel | 0.05 | 1 | AM | EXCL | +0.0448 | 0.0004 | 1.157 |
| Vessel | 0.05 | 1 | AM | Madow | +0.0564 | 0.0004 | 1.198 |
| Vessel | 0.05 | 2 | closed-loop | EXCL | +0.0315 | 0.0004 | 1.130 |
| Vessel | 0.05 | 2 | RS | EXCL | +0.0233 | 0.0004 | 1.110 |
| Vessel | 0.05 | 2 | RS | Madow | +0.0329 | 0.0004 | 1.155 |
| Vessel | 0.05 | 2 | AM | EXCL | +0.0479 | 0.0004 | 1.157 |
| Vessel | 0.05 | 2 | AM | Madow | +0.0609 | 0.0004 | 1.200 |
| Vessel | 0.1 | 0 | closed-loop | EXCL | +0.0570 | 0.0004 | 1.158 |
| Vessel | 0.1 | 0 | RS | EXCL | +0.0449 | 0.0004 | 1.142 |
| Vessel | 0.1 | 0 | RS | Madow | +0.0584 | 0.0004 | 1.185 |
| Vessel | 0.1 | 0 | AM | EXCL | +0.0625 | 0.0004 | 1.143 |
| Vessel | 0.1 | 0 | AM | Madow | +0.0934 | 0.0004 | 1.214 |
| Vessel | 0.1 | 1 | closed-loop | EXCL | +0.0631 | 0.0004 | 1.161 |
| Vessel | 0.1 | 1 | RS | EXCL | +0.0509 | 0.0004 | 1.148 |
| Vessel | 0.1 | 1 | RS | Madow | +0.0679 | 0.0004 | 1.198 |
| Vessel | 0.1 | 1 | AM | EXCL | +0.0656 | 0.0004 | 1.140 |
| Vessel | 0.1 | 1 | AM | Madow | +0.0959 | 0.0004 | 1.204 |
| Vessel | 0.1 | 2 | closed-loop | EXCL | +0.0640 | 0.0004 | 1.155 |
| Vessel | 0.1 | 2 | RS | EXCL | +0.0533 | 0.0004 | 1.146 |
| Vessel | 0.1 | 2 | RS | Madow | +0.0711 | 0.0004 | 1.195 |
| Vessel | 0.1 | 2 | AM | EXCL | +0.0649 | 0.0004 | 1.132 |
| Vessel | 0.1 | 2 | AM | Madow | +0.0948 | 0.0004 | 1.193 |

## 2. Blocks of EXCL

Statistics of the blocks formed by EXCL in one replication per test deployment, pooled over blocks. "Closed by" gives the share of blocks that ended because the next block starts at a recording start, at an onset alignment, or because the mass would exceed 1 (mass cap); the remainder ends with the data. Closed-loop rule and EXCL draws at the probabilities of the two open-loop rules.

| Dataset | b | Rule | Blocks | Units per block (mean) | median | Mass per block (mean) | p5 | p95 | Closed by recording start | Closed by onset alignment | Closed by mass cap | End of data |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | closed-loop | 8504 | 6.20 | 5 | 0.305 | 0.057 | 0.944 | 0.638 | 0.306 | 0.055 | 0.001 |
| Whale | 0.05 | AM+EXCL | 8503 | 6.20 | 5 | 0.304 | 0.042 | 0.901 | 0.638 | 0.306 | 0.055 | 0.001 |
| Whale | 0.05 | RS+EXCL | 8463 | 6.23 | 5 | 0.308 | 0.079 | 0.952 | 0.641 | 0.307 | 0.051 | 0.001 |
| Whale | 0.1 | closed-loop | 9796 | 5.38 | 5 | 0.541 | 0.114 | 0.975 | 0.554 | 0.266 | 0.180 | 0.001 |
| Whale | 0.1 | AM+EXCL | 9879 | 5.34 | 5 | 0.524 | 0.086 | 0.971 | 0.549 | 0.263 | 0.187 | 0.001 |
| Whale | 0.1 | RS+EXCL | 9577 | 5.51 | 5 | 0.544 | 0.141 | 0.973 | 0.567 | 0.272 | 0.161 | 0.001 |
| Vessel | 0.05 | closed-loop | 1758 | 10.68 | 10 | 0.544 | 0.122 | 0.985 | 0.071 | 0.776 | 0.152 | 0.001 |
| Vessel | 0.05 | AM+EXCL | 1804 | 10.41 | 9 | 0.527 | 0.092 | 0.988 | 0.069 | 0.756 | 0.174 | 0.001 |
| Vessel | 0.05 | RS+EXCL | 1744 | 10.77 | 10 | 0.539 | 0.125 | 0.984 | 0.072 | 0.782 | 0.145 | 0.001 |
| Vessel | 0.1 | closed-loop | 2672 | 7.03 | 7 | 0.709 | 0.124 | 0.990 | 0.047 | 0.510 | 0.442 | 0.001 |
| Vessel | 0.1 | AM+EXCL | 2759 | 6.81 | 5 | 0.689 | 0.143 | 0.990 | 0.045 | 0.494 | 0.460 | 0.001 |
| Vessel | 0.1 | RS+EXCL | 2652 | 7.08 | 8 | 0.709 | 0.114 | 0.990 | 0.047 | 0.514 | 0.438 | 0.001 |

## 3. Controller parameters

Seven configurations of the closed-loop rule on test data: the base configuration (p_min = b/4, token-bucket depth 100, step size eta = 0.05) and p_min halved (b/8) or doubled (b/2), depth halved (50) or doubled (200), eta halved (0.025) or doubled (0.10). 2000 replications; cells are means over deployments. Widths are relative to Uniform+PPI of the same run ("nominal": the interval as reported; "matched": scaled to 95% coverage). E[V]/Var is the mean of the variance estimate relative to the empirical variance of the point estimate. The base configuration was run with different random numbers from those of the main tables, so its entries differ from them within Monte-Carlo error.

| Configuration | Dataset | b | Coverage EXCL | Coverage Ind. | Recall Ind. | Recall EXCL | Recall ratio | Rate/b (EXCL) | Nominal width Ind. | Nominal width EXCL | Matched width Ind. | Matched width EXCL | Bias z (EXCL) | E[V]/Var (EXCL) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base | Whale | 0.05 | 0.9600 | 0.9445 | 0.2113 | 0.2357 | 1.116 | 1.018 | 1.092 | 1.113 | 1.112 | 1.046 | -0.35 | 1.183 |
| base | Whale | 0.1 | 0.9723 | 0.9499 | 0.3622 | 0.4317 | 1.192 | 1.008 | 1.101 | 1.149 | 1.107 | 1.006 | +0.30 | 1.330 |
| base | Vessel | 0.05 | 0.9710 | 0.9470 | 0.2258 | 0.2532 | 1.121 | 0.997 | 1.002 | 1.024 | 1.001 | 0.897 | +0.10 | 1.252 |
| base | Vessel | 0.1 | 0.9838 | 0.9523 | 0.3920 | 0.4543 | 1.159 | 0.999 | 1.004 | 1.052 | 1.001 | 0.839 | -0.03 | 1.555 |
| p_min = b/8 | Whale | 0.05 | 0.9639 | 0.9476 | 0.2128 | 0.2377 | 1.117 | 1.019 | 1.108 | 1.129 | 1.132 | 1.068 | +0.21 | 1.214 |
| p_min = b/8 | Whale | 0.1 | 0.9724 | 0.9492 | 0.3645 | 0.4348 | 1.193 | 1.008 | 1.118 | 1.162 | 1.136 | 1.049 | +0.01 | 1.325 |
| p_min = b/8 | Vessel | 0.05 | 0.9748 | 0.9483 | 0.2261 | 0.2534 | 1.121 | 0.998 | 1.003 | 1.024 | 1.013 | 0.889 | +0.82 | 1.300 |
| p_min = b/8 | Vessel | 0.1 | 0.9760 | 0.9495 | 0.3917 | 0.4548 | 1.161 | 0.999 | 1.003 | 1.051 | 1.020 | 0.896 | -0.14 | 1.429 |
| p_min = b/2 | Whale | 0.05 | 0.9675 | 0.9487 | 0.2047 | 0.2277 | 1.113 | 1.014 | 1.059 | 1.081 | 1.051 | 0.982 | -0.32 | 1.241 |
| p_min = b/2 | Whale | 0.1 | 0.9762 | 0.9519 | 0.3527 | 0.4185 | 1.187 | 1.006 | 1.066 | 1.113 | 1.061 | 0.962 | +0.43 | 1.368 |
| p_min = b/2 | Vessel | 0.05 | 0.9712 | 0.9445 | 0.2257 | 0.2525 | 1.119 | 0.998 | 1.003 | 1.023 | 1.034 | 0.910 | -0.43 | 1.273 |
| p_min = b/2 | Vessel | 0.1 | 0.9808 | 0.9457 | 0.3912 | 0.4543 | 1.161 | 0.999 | 1.000 | 1.051 | 1.026 | 0.879 | +1.22 | 1.499 |
| depth 50 | Whale | 0.05 | 0.9635 | 0.9499 | 0.2109 | 0.2349 | 1.114 | 1.017 | 1.108 | 1.124 | 1.086 | 1.022 | -0.28 | 1.197 |
| depth 50 | Whale | 0.1 | 0.9711 | 0.9478 | 0.3628 | 0.4304 | 1.186 | 1.007 | 1.117 | 1.156 | 1.143 | 1.045 | -0.61 | 1.297 |
| depth 50 | Vessel | 0.05 | 0.9652 | 0.9427 | 0.2262 | 0.2531 | 1.119 | 0.997 | 1.004 | 1.024 | 1.016 | 0.922 | -0.65 | 1.251 |
| depth 50 | Vessel | 0.1 | 0.9820 | 0.9495 | 0.3916 | 0.4545 | 1.161 | 0.999 | 1.003 | 1.053 | 1.021 | 0.885 | +0.18 | 1.488 |
| depth 200 | Whale | 0.05 | 0.9640 | 0.9467 | 0.2112 | 0.2358 | 1.116 | 1.018 | 1.091 | 1.112 | 1.090 | 1.028 | +0.26 | 1.206 |
| depth 200 | Whale | 0.1 | 0.9713 | 0.9454 | 0.3621 | 0.4317 | 1.192 | 1.008 | 1.100 | 1.146 | 1.107 | 1.003 | -0.75 | 1.346 |
| depth 200 | Vessel | 0.05 | 0.9692 | 0.9503 | 0.2265 | 0.2533 | 1.118 | 0.998 | 1.007 | 1.025 | 1.036 | 0.938 | +0.46 | 1.252 |
| depth 200 | Vessel | 0.1 | 0.9840 | 0.9432 | 0.3919 | 0.4548 | 1.160 | 0.999 | 1.003 | 1.052 | 1.049 | 0.872 | +0.52 | 1.484 |
| eta 0.025 | Whale | 0.05 | 0.9615 | 0.9491 | 0.2159 | 0.2418 | 1.120 | 1.034 | 1.088 | 1.110 | 1.093 | 1.040 | +0.16 | 1.193 |
| eta 0.025 | Whale | 0.1 | 0.9730 | 0.9488 | 0.3660 | 0.4374 | 1.195 | 1.017 | 1.102 | 1.149 | 1.115 | 1.004 | -0.39 | 1.368 |
| eta 0.025 | Vessel | 0.05 | 0.9720 | 0.9475 | 0.2256 | 0.2526 | 1.120 | 0.994 | 1.002 | 1.023 | 1.008 | 0.921 | -1.21 | 1.237 |
| eta 0.025 | Vessel | 0.1 | 0.9815 | 0.9512 | 0.3910 | 0.4539 | 1.161 | 0.997 | 1.001 | 1.052 | 1.011 | 0.883 | -1.40 | 1.435 |
| eta 0.10 | Whale | 0.05 | 0.9639 | 0.9466 | 0.2084 | 0.2321 | 1.114 | 1.008 | 1.101 | 1.119 | 1.108 | 1.027 | +0.20 | 1.207 |
| eta 0.10 | Whale | 0.1 | 0.9762 | 0.9495 | 0.3602 | 0.4276 | 1.187 | 1.003 | 1.107 | 1.150 | 1.119 | 1.013 | -0.17 | 1.333 |
| eta 0.10 | Vessel | 0.05 | 0.9728 | 0.9512 | 0.2263 | 0.2528 | 1.117 | 0.999 | 1.010 | 1.027 | 0.985 | 0.898 | -1.19 | 1.290 |
| eta 0.10 | Vessel | 0.1 | 0.9798 | 0.9507 | 0.3927 | 0.4545 | 1.157 | 0.999 | 1.009 | 1.052 | 0.999 | 0.867 | -0.43 | 1.442 |

### Block-mass cap

EXCL closes a block when its mass would exceed L (L = 1 in the paper). Closed-loop rule, 2000 replications on test data; cells are means over deployments. The last column counts, per replication, units whose probability exceeds L (such a unit cannot be drawn with its exact probability within a block of mass L).

| Dataset | b | L | Coverage | Lowest deployment coverage | Event recall | Recall SE | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Units per block | Mass per block | Blocks per run | Units with p > L (per run) | Rate/b | Bias z |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | 0.5 | 0.9622 | 0.9525 | 0.2317 | 0.0001 | 1.115 | 1.037 | 5.74 | 0.291 | 1966 | 0.00 | 1.018 | +0.08 |
| Whale | 0.05 | 0.67 | 0.9652 | 0.9580 | 0.2334 | 0.0001 | 1.115 | 1.030 | 6.55 | 0.332 | 1821 | 0.00 | 1.018 | +0.37 |
| Whale | 0.05 | 0.8 | 0.9665 | 0.9555 | 0.2344 | 0.0001 | 1.112 | 1.016 | 7.01 | 0.355 | 1761 | 0.00 | 1.018 | -0.25 |
| Whale | 0.05 | 1 | 0.9653 | 0.9545 | 0.2357 | 0.0001 | 1.113 | 1.025 | 7.48 | 0.379 | 1707 | 0.00 | 1.018 | -0.63 |
| Whale | 0.1 | 0.5 | 0.9667 | 0.9575 | 0.3924 | 0.0001 | 1.148 | 1.050 | 3.44 | 0.346 | 3160 | 17.94 | 1.008 | -0.90 |
| Whale | 0.1 | 0.67 | 0.9722 | 0.9610 | 0.4145 | 0.0001 | 1.149 | 1.030 | 4.45 | 0.448 | 2451 | 0.27 | 1.008 | -0.07 |
| Whale | 0.1 | 0.8 | 0.9706 | 0.9630 | 0.4259 | 0.0001 | 1.150 | 1.024 | 5.10 | 0.514 | 2151 | 0.00 | 1.008 | +0.56 |
| Whale | 0.1 | 1 | 0.9731 | 0.9640 | 0.4318 | 0.0001 | 1.148 | 1.013 | 5.76 | 0.579 | 1964 | 0.00 | 1.008 | +0.14 |
| Vessel | 0.05 | 0.5 | 0.9685 | 0.9675 | 0.2453 | 0.0002 | 1.027 | 0.929 | 7.10 | 0.354 | 1324 | 0.00 | 0.998 | +0.64 |
| Vessel | 0.05 | 0.67 | 0.9680 | 0.9665 | 0.2493 | 0.0002 | 1.025 | 0.927 | 8.73 | 0.435 | 1077 | 0.00 | 0.998 | +0.02 |
| Vessel | 0.05 | 0.8 | 0.9690 | 0.9660 | 0.2516 | 0.0002 | 1.025 | 0.919 | 9.69 | 0.483 | 971 | 0.00 | 0.998 | +0.26 |
| Vessel | 0.05 | 1 | 0.9728 | 0.9690 | 0.2535 | 0.0002 | 1.025 | 0.913 | 10.77 | 0.537 | 873 | 0.00 | 0.998 | +0.36 |
| Vessel | 0.1 | 0.5 | 0.9712 | 0.9680 | 0.4203 | 0.0003 | 1.056 | 0.941 | 3.90 | 0.389 | 2411 | 0.00 | 0.999 | -0.41 |
| Vessel | 0.1 | 0.67 | 0.9728 | 0.9680 | 0.4350 | 0.0003 | 1.054 | 0.913 | 5.11 | 0.511 | 1837 | 0.00 | 0.999 | +0.79 |
| Vessel | 0.1 | 0.8 | 0.9815 | 0.9815 | 0.4439 | 0.0003 | 1.052 | 0.906 | 5.96 | 0.595 | 1578 | 0.00 | 0.999 | -0.31 |
| Vessel | 0.1 | 1 | 0.9825 | 0.9815 | 0.4549 | 0.0003 | 1.053 | 0.875 | 7.10 | 0.709 | 1325 | 0.00 | 0.999 | +1.44 |

## 4. Variance estimates for EXCL

V1 = N^-2 sum_k T_k^2 is the estimate V-hat of Eq. (1) of the paper (EXCL uploads at most one unit per block, so V-hat is the sum of the squared block totals), with T_k = sum over the units t of block k of xi_t r_t / p_t. V2 is the successive-difference estimate over consecutive blocks of one recording, V2 = N^-2 (1/2) sum_k w_k, with w_k = (T_k - T_(k+1))^2 for a block that has a successor in its recording, w_k = (T_(k-1) - T_k)^2 for the last block of a recording with at least two blocks, and w_k = 2 T_k^2 for a recording with a single block. The intervals use max(V, V_plug) as in the paper ("max") or the raw estimate. Columns: coverage of the interval with V1 (paper) and with V2 (both with max), the width ratio V2/V1 (raw and max form), E[V]/Var, and the share of replications in which V_plug exceeds V1 or V2. 2000 replications on test data, EXCL draws at the probabilities of the closed-loop rule, Robust Sampling (RS) and Active Measurement (AM); cells are means over deployments.

| Dataset | b | Rule | Coverage, max(V1, Vp) (paper) | Coverage, max(V2, Vp) | Lowest deployment coverage, max(V2, Vp) | Half-width V2/V1 (raw) | Half-width ratio (max form) | E[V1]/Var | E[V2]/Var | E[V_plug]/Var | Share V_plug > V1 | Share V_plug > V2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | closed-loop | 0.9640 | 0.9604 | 0.9490 | 0.969 | 0.983 | 1.164 | 1.089 | 0.944 | 0.214 | 0.277 |
| Whale | 0.05 | RS+EXCL | 0.9685 | 0.9631 | 0.9525 | 0.964 | 0.976 | 1.225 | 1.135 | 1.016 | 0.133 | 0.227 |
| Whale | 0.05 | AM+EXCL | 0.9577 | 0.9556 | 0.9325 | 0.980 | 0.989 | 1.131 | 1.086 | 0.946 | 0.229 | 0.268 |
| Whale | 0.1 | closed-loop | 0.9715 | 0.9660 | 0.9525 | 0.947 | 0.973 | 1.302 | 1.160 | 1.069 | 0.187 | 0.244 |
| Whale | 0.1 | RS+EXCL | 0.9810 | 0.9758 | 0.9680 | 0.940 | 0.965 | 1.412 | 1.241 | 1.181 | 0.126 | 0.220 |
| Whale | 0.1 | AM+EXCL | 0.9648 | 0.9614 | 0.9485 | 0.961 | 0.982 | 1.210 | 1.111 | 1.033 | 0.211 | 0.250 |
| Vessel | 0.05 | closed-loop | 0.9730 | 0.9695 | 0.9655 | 0.968 | 0.980 | 1.241 | 1.162 | 1.156 | 0.251 | 0.494 |
| Vessel | 0.05 | RS+EXCL | 0.9705 | 0.9657 | 0.9615 | 0.967 | 0.979 | 1.278 | 1.195 | 1.194 | 0.208 | 0.499 |
| Vessel | 0.05 | AM+EXCL | 0.9650 | 0.9623 | 0.9615 | 0.971 | 0.983 | 1.198 | 1.130 | 1.098 | 0.344 | 0.470 |
| Vessel | 0.1 | closed-loop | 0.9818 | 0.9792 | 0.9740 | 0.943 | 0.969 | 1.502 | 1.336 | 1.403 | 0.148 | 0.760 |
| Vessel | 0.1 | RS+EXCL | 0.9848 | 0.9828 | 0.9820 | 0.949 | 0.973 | 1.500 | 1.352 | 1.413 | 0.130 | 0.785 |
| Vessel | 0.1 | AM+EXCL | 0.9688 | 0.9640 | 0.9580 | 0.953 | 0.974 | 1.216 | 1.102 | 1.128 | 0.296 | 0.587 |

## 5. The single-start condition

Condition (ii) of Proposition 2 can be checked on development data. For a fixed predictor, the check quantity is kappa = N^2 Var(Madow) / sum_t r_t^2 / p_t, computed exactly (by scanning the cumulative-probability axis) for each deployment, with the prequential predictor evaluated at its expected weights; by the identity of Proposition 2(ii), N^2 (E V-hat - Var) = (1 - kappa) sum_t r_t^2 / p_t, so the single-start interval is conservative in expectation exactly when kappa <= 1. The table gives the largest and the pooled kappa over the deployments of a cell on development and on test data, together with the coverage of the single-start interval on test data (taken from the Madow runs of the main tables). RS = Robust Sampling, AM = Active Measurement.

| Dataset | b | Rule | Development deployments | kappa dev, largest | Development deployments with kappa > 1 | kappa dev, pooled | kappa test, largest | kappa test, pooled | Single-start coverage, test (mean) | Single-start coverage, test (lowest deployment) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale | 0.05 | RS | 4 | 0.760 | 0 | 0.716 | 0.875 | 0.788 | 0.9749 | 0.9660 |
| Whale | 0.05 | AM | 4 | 0.847 | 0 | 0.815 | 0.936 | 0.872 | 0.9643 | 0.9510 |
| Whale | 0.1 | RS | 4 | 0.714 | 0 | 0.660 | 0.750 | 0.685 | 0.9846 | 0.9750 |
| Whale | 0.1 | AM | 4 | 0.820 | 0 | 0.777 | 0.874 | 0.815 | 0.9693 | 0.9555 |
| Vessel | 0.05 | RS | 2 | 0.899 | 0 | 0.838 | 0.712 | 0.711 | 0.9772 | 0.9740 |
| Vessel | 0.05 | AM | 2 | 0.809 | 0 | 0.803 | 0.838 | 0.819 | 0.9725 | 0.9720 |
| Vessel | 0.1 | RS | 2 | 0.762 | 0 | 0.700 | 0.634 | 0.623 | 0.9880 | 0.9855 |
| Vessel | 0.1 | AM | 2 | 0.735 | 0 | 0.697 | 0.773 | 0.708 | 0.9770 | 0.9680 |
| Gray's Reef | 0.05 | RS | 2 | 1.094 | 1 | 0.824 | 0.987 | 0.847 | 0.9686 | 0.9465 |
| Gray's Reef | 0.05 | AM | 2 | 1.223 | 1 | 0.903 | 1.131 | 0.804 | 0.9841 | 0.9750 |
| Gray's Reef | 0.1 | RS | 2 | 1.003 | 1 | 0.963 | 0.689 | 0.499 | 0.9921 | 0.9850 |
| Gray's Reef | 0.1 | AM | 2 | 0.860 | 0 | 0.826 | 0.962 | 0.773 | 0.9853 | 0.9790 |

Gray's Reef development deployments with kappa > 1: 3 of 8 (2 deployments x 2 rules x 2 budgets). On the whale and vessel data kappa <= 1 on every development and test deployment (largest value on test data: 0.936).

### Periodic stress test

A residual component a cos(2 pi C_t) is added to the prequential predictor, where C_t is the within-recording cumulative probability (so the residual is periodic with the sampling interval of the single start). Vessel test data, Active Measurement probabilities, b = 0.05, 2000 replications; cells are means over the two stations. Madow and EXCL use the same probabilities; kappa is the exact check quantity, shared by both draws.

| Draw | Amplitude a | kappa (largest) | kappa (pooled) | Coverage | Lowest station coverage | Bias z | Matched width (absolute) | E[V-hat]/Var |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Madow (single start) | 0 | 0.838 | 0.819 | 0.9705 | 0.9670 | +0.64 | 0.0730 | 1.228 |
| Madow (single start) | 0.05 | 0.868 | 0.863 | 0.9672 | 0.9665 | +0.50 | 0.0754 | 1.167 |
| Madow (single start) | 0.1 | 1.015 | 0.995 | 0.9535 | 0.9525 | +0.40 | 0.0809 | 1.018 |
| Madow (single start) | 0.2 | 1.443 | 1.397 | 0.9070 | 0.8980 | +0.35 | 0.1005 | 0.733 |
| EXCL | 0 | 0.838 | 0.819 | 0.9645 | 0.9625 | -0.92 | 0.0758 | 1.155 |
| EXCL | 0.05 | 0.868 | 0.863 | 0.9643 | 0.9640 | -0.80 | 0.0759 | 1.151 |
| EXCL | 0.1 | 1.015 | 0.995 | 0.9677 | 0.9655 | -0.71 | 0.0762 | 1.148 |
| EXCL | 0.2 | 1.443 | 1.397 | 0.9698 | 0.9650 | -0.59 | 0.0791 | 1.147 |
