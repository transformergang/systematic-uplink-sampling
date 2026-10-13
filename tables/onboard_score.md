# Vessel data with the integer on-board score

The vessel result is repeated with the score computed by the integer chain that runs on the node (25--500 Hz band-pass at 2 kHz, log, causal baseline, binning; see `code/node/score/`) instead of the 14-band score used for the main tables. The chain was evaluated offline on the stored audio; clips shorter than 20000 decimated samples receive no score and are dropped (2.1% of the clips), so the unit sets of the two scores differ marginally. The price level kappa is re-selected on development data for each score (0.80 for the 14-band score, 0.90 for the integer chain). Closed-loop rule, 2000 replications on the two test stations; cells are means over the two stations, recall is pooled over events.

| Score | AUC (test) | b | Coverage EXCL | Coverage Ind. | Recall EXCL | Recall Ind. | Recall ratio | Matched width EXCL / Uniform+PPI | Matched width Ind. / Uniform+PPI | Matched width EXCL / Ind. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 14-band score | 0.769 | 0.05 | 0.971 | 0.948 | 0.253 | 0.226 | 1.119 | 0.945 | 1.022 | 0.924 |
| 14-band score | 0.769 | 0.10 | 0.977 | 0.952 | 0.454 | 0.392 | 1.159 | 0.880 | 0.995 | 0.885 |
| integer on-board score | 0.710 | 0.05 | 0.963 | 0.950 | 0.254 | 0.227 | 1.118 | 0.979 | 1.016 | 0.964 |
| integer on-board score | 0.710 | 0.10 | 0.983 | 0.949 | 0.452 | 0.394 | 1.148 | 0.921 | 1.058 | 0.870 |
