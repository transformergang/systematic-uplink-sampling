# Main results on test data

Every cell of the paper's main comparison. Sampling is replicated 2000 times on each fixed test population; a cell is the mean over the test deployments of a dataset (coverage, matched width, upload rate) with event recall pooled over events. *Coverage* is that of the nominal 95% interval; *matched width* is the half-width scaled by the smallest common factor that gives 95% coverage over the replications, relative to Uniform+PPI (an oracle diagnostic of actual error); *event recall* is the fraction of events with at least one upload; *upload rate / b* is the realised fraction of uploaded units divided by the budget. The union bound is the largest event recall any design with these marginals can reach. Trigger was run on the whale and vessel data only. *Nominal width* is the mean width of the interval as reported, relative to Uniform+PPI. The price parameter kappa of the closed-loop rule, selected on development data, is 0.95 (whale), 0.80 (vessel) and 0 (Gray's Reef).

## Closed-loop rule

Designs: Uniform+PPI (probability b), Duty cycle, Independent (the same closed-loop rule with independent draws and per-unit price updates), EXCL (blocked restart), EXCL without onset alignment (blocks start only at a recording start or when the mass would exceed 1), Madow with a single start per recording, where the price is driven either by the block mass ("open-loop price", exact marginals) or by realised uploads ("closed-loop price", the marginals are no longer exact).

### Whale calls (BioDCASE 2026 Task 2), b = 0.05

| Design | Coverage | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Event recall | Upload rate / b |
|---|---:|---:|---:|---:|---:|
| Trigger (top scores) | 0.000 | -- | -- | 0.226 | 1.059 |
| Uniform+PPI | 0.953 | 1.000 | 1.000 | 0.171 | 0.999 |
| Duty cycle | 0.980 | 1.056 | 0.787 | 0.189 | 1.000 |
| Independent, same rule | 0.945 | 1.093 | 1.128 | 0.211 | 1.018 |
| EXCL (blocked restart) | 0.962 | 1.113 | 1.048 | 0.236 | 1.018 |
| EXCL without onset alignment | 0.969 | 1.112 | 1.020 | 0.236 | 1.018 |
| Madow, open-loop price | 0.970 | 1.113 | 1.008 | 0.239 | 1.018 |
| Madow, closed-loop price | 0.971 | 1.112 | 1.005 | 0.238 | 1.018 |
| Union bound, same $p_t$ | -- | -- | -- | 0.239 | -- |


### Whale calls (BioDCASE 2026 Task 2), b = 0.10

| Design | Coverage | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Event recall | Upload rate / b |
|---|---:|---:|---:|---:|---:|
| Trigger (top scores) | 0.000 | -- | -- | 0.381 | 1.054 |
| Uniform+PPI | 0.949 | 1.000 | 1.000 | 0.305 | 0.999 |
| Duty cycle | 1.000 | 0.997 | 0.766 | 0.357 | 1.000 |
| Independent, same rule | 0.952 | 1.101 | 1.081 | 0.362 | 1.008 |
| EXCL (blocked restart) | 0.973 | 1.148 | 1.003 | 0.432 | 1.008 |
| EXCL without onset alignment | 0.976 | 1.147 | 0.978 | 0.432 | 1.008 |
| Madow, open-loop price | 0.975 | 1.146 | 0.970 | 0.437 | 1.008 |
| Madow, closed-loop price | 0.977 | 1.146 | 0.968 | 0.437 | 1.008 |
| Union bound, same $p_t$ | -- | -- | -- | 0.438 | -- |


### Vessel presence (VLIZ), b = 0.05

| Design | Coverage | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Event recall | Upload rate / b |
|---|---:|---:|---:|---:|---:|
| Trigger (top scores) | 0.000 | -- | -- | 0.372 | 1.133 |
| Uniform+PPI | 0.953 | 1.000 | 1.000 | 0.178 | 1.000 |
| Duty cycle | 0.950 | 1.083 | 0.953 | 0.201 | 1.000 |
| Independent, same rule | 0.948 | 1.004 | 1.022 | 0.226 | 0.998 |
| EXCL (blocked restart) | 0.971 | 1.025 | 0.945 | 0.253 | 0.998 |
| EXCL without onset alignment | 0.973 | 1.024 | 0.926 | 0.252 | 0.997 |
| Madow, open-loop price | 0.980 | 1.022 | 0.887 | 0.261 | 0.997 |
| Madow, closed-loop price | 0.980 | 1.027 | 0.885 | 0.262 | 0.997 |
| Union bound, same $p_t$ | -- | -- | -- | 0.261 | -- |


### Vessel presence (VLIZ), b = 0.10

| Design | Coverage | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Event recall | Upload rate / b |
|---|---:|---:|---:|---:|---:|
| Trigger (top scores) | 0.000 | -- | -- | 0.535 | 1.114 |
| Uniform+PPI | 0.951 | 1.000 | 1.000 | 0.318 | 1.000 |
| Duty cycle | 1.000 | 1.075 | 0.652 | 0.379 | 1.000 |
| Independent, same rule | 0.952 | 1.004 | 0.995 | 0.392 | 0.999 |
| EXCL (blocked restart) | 0.977 | 1.051 | 0.880 | 0.454 | 0.999 |
| EXCL without onset alignment | 0.986 | 1.051 | 0.844 | 0.448 | 0.999 |
| Madow, open-loop price | 0.988 | 1.050 | 0.814 | 0.475 | 0.999 |
| Madow, closed-loop price | 0.983 | 1.046 | 0.828 | 0.473 | 0.999 |
| Union bound, same $p_t$ | -- | -- | -- | 0.479 | -- |


### Gray's Reef (SanctSound), b = 0.05

| Design | Coverage | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Event recall | Upload rate / b |
|---|---:|---:|---:|---:|---:|
| Trigger (top scores) | not run | -- | -- | -- | -- |
| Uniform+PPI | 0.942 | 1.000 | 1.000 | 0.079 | 1.000 |
| Duty cycle | 0.937 | 0.928 | 0.875 | 0.083 | 1.000 |
| Independent, same rule | 0.949 | 0.908 | 0.871 | 0.124 | 1.008 |
| EXCL (blocked restart) | 0.969 | 0.932 | 0.805 | 0.131 | 1.008 |
| EXCL without onset alignment | 0.970 | 0.933 | 0.779 | 0.131 | 1.008 |
| Madow, open-loop price | 0.980 | 0.930 | 0.709 | 0.133 | 1.008 |
| Madow, closed-loop price | 1.000 | 0.924 | 0.534 | 0.137 | 1.007 |
| Union bound, same $p_t$ | -- | -- | -- | 0.133 | -- |


### Gray's Reef (SanctSound), b = 0.10

| Design | Coverage | Nominal width / Uniform+PPI | Matched width / Uniform+PPI | Event recall | Upload rate / b |
|---|---:|---:|---:|---:|---:|
| Trigger (top scores) | not run | -- | -- | -- | -- |
| Uniform+PPI | 0.955 | 1.000 | 1.000 | 0.152 | 1.000 |
| Duty cycle | 0.975 | 0.923 | 0.735 | 0.163 | 1.000 |
| Independent, same rule | 0.953 | 0.907 | 0.916 | 0.232 | 1.004 |
| EXCL (blocked restart) | 0.980 | 0.961 | 0.822 | 0.246 | 1.004 |
| EXCL without onset alignment | 0.981 | 0.961 | 0.797 | 0.245 | 1.004 |
| Madow, open-loop price | 0.991 | 0.958 | 0.700 | 0.251 | 1.004 |
| Madow, closed-loop price | 1.000 | 0.954 | 0.539 | 0.246 | 1.004 |
| Union bound, same $p_t$ | -- | -- | -- | 0.253 | -- |


## Open-loop rules (Robust Sampling, Active Measurement) with three draws

The probabilities are those of the open-loop rule; only the draw changes. "Union bound" is the event-recall bound at these probabilities (the same for the three draws).

### Whale calls (BioDCASE 2026 Task 2), b = 0.05

| Rule | Draw | Coverage | Matched width / Uniform+PPI | Event recall | Upload rate / b | Union bound |
|---|---:|---:|---:|---:|---:|---:|
| Robust Sampling | independent | 0.952 | 0.945 | 0.185 | 0.993 | 0.207 |
| Robust Sampling | EXCL | 0.970 | 0.881 | 0.204 | 0.993 |  |
| Robust Sampling | Madow (single start) | 0.975 | 0.857 | 0.206 | 0.993 |  |
| Active Measurement | independent | 0.951 | 1.318 | 0.228 | 0.972 | 0.263 |
| Active Measurement | EXCL | 0.958 | 1.291 | 0.260 | 0.971 |  |
| Active Measurement | Madow (single start) | 0.964 | 1.232 | 0.262 | 0.971 |  |


### Whale calls (BioDCASE 2026 Task 2), b = 0.10

| Rule | Draw | Coverage | Matched width / Uniform+PPI | Event recall | Upload rate / b | Union bound |
|---|---:|---:|---:|---:|---:|---:|
| Robust Sampling | independent | 0.949 | 0.946 | 0.327 | 0.994 | 0.388 |
| Robust Sampling | EXCL | 0.980 | 0.834 | 0.382 | 0.993 |  |
| Robust Sampling | Madow (single start) | 0.985 | 0.826 | 0.386 | 0.993 |  |
| Active Measurement | independent | 0.951 | 1.319 | 0.388 | 0.971 | 0.482 |
| Active Measurement | EXCL | 0.965 | 1.250 | 0.472 | 0.971 |  |
| Active Measurement | Madow (single start) | 0.969 | 1.226 | 0.481 | 0.971 |  |


### Vessel presence (VLIZ), b = 0.05

| Rule | Draw | Coverage | Matched width / Uniform+PPI | Event recall | Upload rate / b | Union bound |
|---|---:|---:|---:|---:|---:|---:|
| Robust Sampling | independent | 0.953 | 0.958 | 0.197 | 1.001 | 0.224 |
| Robust Sampling | EXCL | 0.970 | 0.898 | 0.217 | 1.001 |  |
| Robust Sampling | Madow (single start) | 0.977 | 0.857 | 0.224 | 1.001 |  |
| Active Measurement | independent | 0.952 | 1.173 | 0.285 | 1.010 | 0.343 |
| Active Measurement | EXCL | 0.960 | 1.137 | 0.331 | 1.011 |  |
| Active Measurement | Madow (single start) | 0.972 | 1.040 | 0.342 | 1.010 |  |


### Vessel presence (VLIZ), b = 0.10

| Rule | Draw | Coverage | Matched width / Uniform+PPI | Event recall | Upload rate / b | Union bound |
|---|---:|---:|---:|---:|---:|---:|
| Robust Sampling | independent | 0.951 | 0.982 | 0.343 | 1.001 | 0.414 |
| Robust Sampling | EXCL | 0.981 | 0.862 | 0.395 | 1.002 |  |
| Robust Sampling | Madow (single start) | 0.988 | 0.825 | 0.412 | 1.001 |  |
| Active Measurement | independent | 0.946 | 1.229 | 0.469 | 1.009 | 0.569 |
| Active Measurement | EXCL | 0.976 | 1.089 | 0.536 | 1.010 |  |
| Active Measurement | Madow (single start) | 0.977 | 1.047 | 0.566 | 1.010 |  |


### Gray's Reef (SanctSound), b = 0.05

| Rule | Draw | Coverage | Matched width / Uniform+PPI | Event recall | Upload rate / b | Union bound |
|---|---:|---:|---:|---:|---:|---:|
| Robust Sampling | independent | 0.945 | 0.984 | 0.082 | 1.003 | 0.087 |
| Robust Sampling | EXCL | 0.963 | 0.910 | 0.086 | 1.002 |  |
| Robust Sampling | Madow (single start) | 0.969 | 0.888 | 0.087 | 1.003 |  |
| Active Measurement | independent | 0.955 | 0.893 | 0.187 | 1.087 | 0.205 |
| Active Measurement | EXCL | 0.971 | 0.839 | 0.198 | 1.087 |  |
| Active Measurement | Madow (single start) | 0.984 | 0.767 | 0.202 | 1.087 |  |


### Gray's Reef (SanctSound), b = 0.10

| Rule | Draw | Coverage | Matched width / Uniform+PPI | Event recall | Upload rate / b | Union bound |
|---|---:|---:|---:|---:|---:|---:|
| Robust Sampling | independent | 0.943 | 0.986 | 0.159 | 1.003 | 0.171 |
| Robust Sampling | EXCL | 0.979 | 0.829 | 0.169 | 1.003 |  |
| Robust Sampling | Madow (single start) | 0.992 | 0.689 | 0.171 | 1.003 |  |
| Active Measurement | independent | 0.955 | 0.903 | 0.338 | 1.087 | 0.367 |
| Active Measurement | EXCL | 0.979 | 0.818 | 0.351 | 1.087 |  |
| Active Measurement | Madow (single start) | 0.985 | 0.798 | 0.365 | 1.087 |  |


## Bias of the point estimate

Standardised bias z = (mean estimate - truth) / (standard error of the mean over the 2000 replications), averaged over the test deployments of a dataset. EXCL and the independent draws are unbiased for every cell; the single start with a closed-loop price is not (Gray's Reef), as Proposition 2(iii) predicts. "Madow" in the open-loop columns uses the open-loop rule's probabilities.

| Dataset | b | Indep. | EXCL | EXCL w/o onset | Madow open-loop price | Madow closed-loop price | RS+EXCL | RS+Madow | AM+EXCL | AM+Madow |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale calls | 0.05 | 0.88 | 0.04 | 0.29 | -0.05 | 0.40 | 0.30 | 1.02 | -0.09 | 0.29 |
| Whale calls | 0.10 | -0.22 | -0.43 | -0.08 | -0.38 | -0.87 | 0.27 | -0.21 | -0.34 | -0.70 |
| Vessel presence | 0.05 | 0.41 | -0.17 | 0.24 | -0.30 | 4.96 | -0.22 | 0.25 | 1.04 | -0.00 |
| Vessel presence | 0.10 | 0.53 | -0.23 | 0.54 | 0.18 | -1.56 | 0.81 | -0.06 | 0.35 | 0.08 |
| Gray's Reef | 0.05 | -0.01 | -0.48 | -0.39 | -0.50 | 10.80 | 0.73 | 0.24 | -0.51 | -0.65 |
| Gray's Reef | 0.10 | 0.57 | -0.27 | -0.19 | 0.35 | -22.69 | 0.17 | 1.68 | -0.03 | -0.39 |

## Trigger: bias of the sample mean

Relative bias of the mean of the uploaded labels (top-scoring segments) as an estimate of the prevalence: (mean estimate - truth) / truth.

| Dataset | b | Coverage | Relative bias | Upload rate / b |
|---|---:|---:|---:|---:|
| Whale calls | 0.05 | 0.000 | +161 % | 1.059 |
| Whale calls | 0.10 | 0.000 | +132 % | 1.054 |
| Vessel presence | 0.05 | 0.000 | +201 % | 1.133 |
| Vessel presence | 0.10 | 0.000 | +165 % | 1.114 |

## Test deployments

Population summary of the test data used in all tables: number of units N, prevalence theta, number of events (maximal runs of positive units, gaps of one negative unit merged; for Gray's Reef the events of the energy detector) and the median event length in units.

| Dataset | Deployment | N | theta | Events | Median event length |
|---|---:|---:|---:|---:|---:|
| Whale | casey2014 | 11631 | 0.422 | 911 | 2.0 |
| Whale | elephantisland2013 | 11188 | 0.524 | 1540 | 5.0 |
| Whale | elephantisland2014 | 12921 | 0.642 | 2143 | 5.0 |
| Whale | kerguelen2005 | 12000 | 0.138 | 599 | 1.0 |
| Whale | maudrise2014 | 5000 | 0.326 | 427 | 2.0 |
| Vessel | vliz_GardenCity | 8741 | 0.182 | 448 | 3.0 |
| Vessel | vliz_Grafton | 10043 | 0.298 | 661 | 3.0 |
| Gray's Reef | sanctsound_gr01_01 | 19869 | 0.007 | 90 | 1.0 |
| Gray's Reef | sanctsound_gr01_02 | 18081 | 0.019 | 192 | 1.0 |
| Gray's Reef | sanctsound_gr01_05 | 20026 | 0.008 | 108 | 1.0 |
| Gray's Reef | sanctsound_gr01_06 | 16692 | 0.018 | 170 | 1.0 |

## Ordered pivotal (Deville) design in the closed loop

The ordered pivotal design (Deville) draws at the probabilities of the closed-loop rule (price driven by the realised uploads) and fixes the fate of a unit only after later units have been seen. Standardised bias z of the estimate for each test deployment and budget (2000 replications), for the ordered pivotal design and for EXCL (the run of the main tables), and the delay between a unit and its decision.

| Deployment / b | z, ordered pivotal | z, EXCL |
|---|---:|---:|
| whale / casey2014 / 0.05 | -0.33 | 0.39 |
| whale / casey2014 / 0.1 | -0.94 | -0.80 |
| whale / elephantisland2013 / 0.05 | -0.61 | 0.18 |
| whale / elephantisland2013 / 0.1 | +1.18 | -0.20 |
| whale / elephantisland2014 / 0.05 | -1.04 | -0.20 |
| whale / elephantisland2014 / 0.1 | -0.04 | -1.46 |
| whale / kerguelen2005 / 0.05 | -1.55 | 0.54 |
| whale / kerguelen2005 / 0.1 | +0.54 | 0.27 |
| whale / maudrise2014 / 0.05 | +1.16 | -0.74 |
| whale / maudrise2014 / 0.1 | +2.00 | 0.06 |
| vliz / vliz_GardenCity / 0.05 | -0.93 | -1.90 |
| vliz / vliz_GardenCity / 0.1 | +1.91 | -0.28 |
| vliz / vliz_Grafton / 0.05 | -0.88 | 1.56 |
| vliz / vliz_Grafton / 0.1 | -0.48 | -0.18 |
| gr / sanctsound_gr01_01 / 0.05 | +1.49 | 0.16 |
| gr / sanctsound_gr01_01 / 0.1 | -0.62 | 0.81 |
| gr / sanctsound_gr01_02 / 0.05 | -0.89 | -1.39 |
| gr / sanctsound_gr01_02 / 0.1 | -1.34 | -0.08 |
| gr / sanctsound_gr01_05 / 0.05 | -1.33 | -0.21 |
| gr / sanctsound_gr01_05 / 0.1 | +1.42 | -0.36 |
| gr / sanctsound_gr01_06 / 0.05 | +1.07 | -0.50 |
| gr / sanctsound_gr01_06 / 0.1 | -0.32 | -1.44 |

Largest |z| over the 22 deployment-budget cells: 2.00 for the ordered pivotal design, 1.90 for EXCL.

| Dataset | Design | Median delay (units) | 99th percentile | Maximum | Mean |
|---|---:|---:|---:|---:|---:|
| Whale | ordered pivotal | 0 | 11 | 60 | 1.000 |
| Whale | ordered pivotal, open-loop price | 0 | 11 | 60 | 1.000 |
| Vessel | ordered pivotal | 0 | 14 | 59 | 1.000 |
| Vessel | ordered pivotal, open-loop price | 0 | 14 | 61 | 1.000 |
| Gray's Reef | ordered pivotal | 0 | 14 | 71 | 1.000 |
| Gray's Reef | ordered pivotal, open-loop price | 0 | 14 | 78 | 1.000 |
