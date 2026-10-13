# Dependence on score quality

The whole procedure (bins, bin statistics, price level kappa chosen on development data, then test data) is re-run on seven versions of the score of each dataset. Let r be the original score converted to its empirical distribution function on the development data (values in [0,1]), u an independent uniform random number per unit, and y the truth of the unit. The levels are N100 (u), N75 (0.25 r + 0.75 u), N50 (0.5 r + 0.5 u), N25 (0.75 r + 0.25 u), L0 (r, the original score), O15 (0.85 r + 0.15 y) and O30 (0.70 r + 0.30 y). O15 and O30 use the truth of the same units and are synthetic, better-than-deployable scorers. AUC is the area under the ROC curve of the score on the test data. Closed-loop rule, 2000 replications; cells are means over test deployments, recall is pooled over events. "Recall ratio" is event recall of EXCL divided by that of independent draws; "matched width EXCL/Ind." is the ratio of matched widths of the two designs (smaller than 1 favours EXCL).

| Dataset | Level | b | kappa | AUC (test) | Coverage EXCL | Coverage Ind. | Recall EXCL | Recall Ind. | Recall ratio | Matched width EXCL/Ind. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Whale | N100 (pure noise) | 0.05 | 0.95 | 0.496 | 0.973 | 0.951 | 0.183 | 0.171 | 1.072 | 0.921 |
| Whale | N100 (pure noise) | 0.10 | 0.95 | 0.496 | 0.982 | 0.955 | 0.345 | 0.306 | 1.129 | 0.865 |
| Whale | N75 (0.25 score + 0.75 noise) | 0.05 | 0.90 | 0.581 | 0.967 | 0.946 | 0.197 | 0.183 | 1.078 | 0.920 |
| Whale | N75 (0.25 score + 0.75 noise) | 0.10 | 0.90 | 0.581 | 0.980 | 0.948 | 0.368 | 0.323 | 1.140 | 0.867 |
| Whale | N50 (0.5 score + 0.5 noise) | 0.05 | 0.90 | 0.700 | 0.966 | 0.948 | 0.216 | 0.196 | 1.103 | 0.922 |
| Whale | N50 (0.5 score + 0.5 noise) | 0.10 | 0.90 | 0.700 | 0.975 | 0.946 | 0.401 | 0.342 | 1.173 | 0.881 |
| Whale | N25 (0.75 score + 0.25 noise) | 0.05 | 0.95 | 0.803 | 0.961 | 0.948 | 0.238 | 0.213 | 1.116 | 0.952 |
| Whale | N25 (0.75 score + 0.25 noise) | 0.10 | 0.95 | 0.803 | 0.974 | 0.954 | 0.434 | 0.365 | 1.190 | 0.940 |
| Whale | L0 (original score) | 0.05 | 0.95 | 0.840 | 0.962 | 0.945 | 0.236 | 0.211 | 1.116 | 0.935 |
| Whale | L0 (original score) | 0.10 | 0.95 | 0.840 | 0.975 | 0.952 | 0.432 | 0.362 | 1.192 | 0.916 |
| Whale | O15 (0.85 score + 0.15 truth) | 0.05 | 0.95 | 0.952 | 0.967 | 0.954 | 0.219 | 0.198 | 1.106 | 0.949 |
| Whale | O15 (0.85 score + 0.15 truth) | 0.10 | 0.95 | 0.952 | 0.976 | 0.957 | 0.395 | 0.342 | 1.153 | 0.930 |
| Whale | O30 (0.70 score + 0.30 truth) | 0.05 | 0.95 | 0.992 | 0.976 | 0.952 | 0.156 | 0.146 | 1.075 | 0.902 |
| Whale | O30 (0.70 score + 0.30 truth) | 0.10 | 0.95 | 0.992 | 0.985 | 0.952 | 0.267 | 0.256 | 1.043 | 0.912 |
| Vessel | N100 (pure noise) | 0.05 | 0.95 | 0.499 | 0.970 | 0.947 | 0.193 | 0.179 | 1.076 | 0.900 |
| Vessel | N100 (pure noise) | 0.10 | 0.95 | 0.499 | 0.977 | 0.948 | 0.356 | 0.319 | 1.115 | 0.888 |
| Vessel | N75 (0.25 score + 0.75 noise) | 0.05 | 0.90 | 0.576 | 0.969 | 0.954 | 0.217 | 0.197 | 1.101 | 0.945 |
| Vessel | N75 (0.25 score + 0.75 noise) | 0.10 | 0.90 | 0.576 | 0.981 | 0.948 | 0.395 | 0.348 | 1.135 | 0.856 |
| Vessel | N50 (0.5 score + 0.5 noise) | 0.05 | 0.80 | 0.682 | 0.964 | 0.948 | 0.235 | 0.213 | 1.105 | 0.942 |
| Vessel | N50 (0.5 score + 0.5 noise) | 0.10 | 0.80 | 0.682 | 0.978 | 0.949 | 0.425 | 0.372 | 1.143 | 0.893 |
| Vessel | N25 (0.75 score + 0.25 noise) | 0.05 | 0.80 | 0.754 | 0.969 | 0.944 | 0.259 | 0.231 | 1.118 | 0.892 |
| Vessel | N25 (0.75 score + 0.25 noise) | 0.10 | 0.80 | 0.754 | 0.980 | 0.946 | 0.462 | 0.399 | 1.157 | 0.862 |
| Vessel | L0 (original score) | 0.05 | 0.80 | 0.769 | 0.973 | 0.956 | 0.253 | 0.226 | 1.119 | 0.953 |
| Vessel | L0 (original score) | 0.10 | 0.80 | 0.769 | 0.983 | 0.953 | 0.455 | 0.392 | 1.162 | 0.863 |
| Vessel | O15 (0.85 score + 0.15 truth) | 0.05 | 0.95 | 0.881 | 0.976 | 0.953 | 0.244 | 0.219 | 1.114 | 0.892 |
| Vessel | O15 (0.85 score + 0.15 truth) | 0.10 | 0.95 | 0.881 | 0.987 | 0.954 | 0.431 | 0.377 | 1.143 | 0.862 |
| Vessel | O30 (0.70 score + 0.30 truth) | 0.05 | 0.95 | 0.964 | 0.981 | 0.951 | 0.254 | 0.229 | 1.111 | 0.844 |
| Vessel | O30 (0.70 score + 0.30 truth) | 0.10 | 0.95 | 0.964 | 0.991 | 0.955 | 0.428 | 0.388 | 1.104 | 0.825 |
