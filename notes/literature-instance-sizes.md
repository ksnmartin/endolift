# Instance sizes in the decision-dependent-probability facility location literature (checked 2026-10-31)

Sources read from the PDFs (pdftotext), not from abstracts.

| paper | model of decision dependence | method | largest solved to optimality | notes |
|---|---|---|---|---|
| Pantuso, arXiv:2509.09816 / J. Comb. Optim. 2025 | finite set of distributions indexed by which *zones* have an open facility (2^|Z| distributions); truncated-normal demand, |S_d| = 50 or 100 scenarios per distribution | decision-dependent L-shaped (Pantuso & Hewitt), CPLEX 20.1 B&B with cuts at integer nodes, single thread, 1800 s limit | all 1,008 instances per demand type with |I| in {10,15,20,25}, |J| in {50,100}, |Z| in {5,7,10}; avg time < 4 min, < 2,000 nodes | stress test |I|=20, |J|=50, |S|=50, 36,000 s limit: |Z|=12 solved (431 s avg), |Z|=15 solved (4,258 s avg), |Z|=17 74/84 within 0.5%, |Z|=20 none (Table 5). So the limit is about 15-17 *influencing* zones at 20 candidate sites. |
| Basciftci, Ahmed, Shen, EJOR 2021 (arXiv:1912.05577) | moment-based DRO: mean and variance of demand piecewise-linear in y; products linearized by McCormick with dual-variable bounds 100 | monolithic MILP in Gurobi 7.5.2 | |I| in 5..10, |J| = 2|I|; 111 s avg at 10x20 (Table 7) | no instances above 10 facilities reported |
| Hellemo, Barton, Tomasgard, CMS 2018 | p_s(y) via exact expression or rational approximation | BARON, selective branching | small illustrative instances (benchmarks Table 2; not read in full, Springer paywalled) | |
| Lugani & Pantuso, Zenodo 10.5281/zenodo.23594350 (Oct 2026) | facility location + capacity planning, endogenous demand; B&B | 5 runs' data only, no paper linked yet | | |

Reading for our paper: the comparable exact-method result is Pantuso 2025: 20-25 candidate sites with 10-15 zones
(2^10 - 2^15 distributions) in minutes to two hours, 17 zones partly, 20 zones not at all in 10 h. Our dense
instances have every one of n = 25-50 sites influencing demand (affine, so "2^n distributions" in his
terminology) and are proven in 17 s to 37 min; the sparse 500-site instances have about 10 influencing sites,
which is inside his solvable zone count but at 20-25x the number of candidate sites. The models differ (affine
p_s(y) over a fixed scenario set vs. a distribution per subset of zones with its own scenario set), so the
comparison is of scale and structure, not head-to-head.
