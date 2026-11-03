# Variational principles (statistical physics / information theory / variational inference) for the decision-dependent expectation sum_s p_s(y) Q_s(y)

Standing notation used throughout (fixed for all sections):

- `y ∈ {0,1}^n` facility-opening vector; `s = 1..S` finite demand scenarios; `Δ` = probability simplex on scenarios.
- `Q_s(y) = min { c^T x : W x ≥ h_s − T_s y, x ≥ 0 }`, convex piecewise-linear in `y`, assumed finite (relatively complete recourse).
- Gibbs/Boltzmann model of decision-dependent probabilities: `p_s(y) = q_s exp(−β E_s(y)) / Z(y)`, `Z(y) = Σ_t q_t exp(−β E_t(y))`, with reference distribution `q ∈ Δ` and "energy" `E_s(y)`. The *affine-energy* case is `E_s(y) = a_s^T y + b_s`.
- Objective: `min_y c_1^T y + Σ_s p_s(y) Q_s(y)`.

The notes label every statement as **EXACT**, **BOUND**, or **APPROXIMATION**, and say whether convexity in `y` survives. Derivations marked "derived here" are my own and were not found stated in the literature; they are elementary consequences of the cited identities but should be checked before use.

---

## KQ1. Gibbs variational principle / Donsker–Varadhan: can Σ_s p_s(y) Q_s(y) be written as a min–max or min–min over a fixed-measure variable r plus entropy, with y entering linearly?

### Takeaway
Donsker–Varadhan represents the **log-partition / entropic-risk functional** as a sup over a fixed simplex, not the expectation itself; the expectation is the *derivative* of that functional, so it has no direct sup/inf representation. However, the DV/Gibbs optimality gap `KL(r‖p(y)) = KL(r‖q) + β E_r[E(y)] + log Z(y)` is a jointly convex function of `(r, y)` once the bilinear term `r_s y_j` is linearized, which for *binary* `y` is exact (McCormick). This turns the decision-dependent two-stage program with affine-energy Gibbs probabilities into an **exact mixed-integer convex (exponential-cone) program** whose continuous relaxation is convex (derived here).

### Cited Findings
- Donsker–Varadhan / Gibbs variational formula, exact statement: for a probability measure `P` and bounded measurable `g`, `log E_P[exp(g)] = sup_{Q ≪ P} { E_Q[g] − KL(Q‖P) }`, equivalently `KL(Q‖P) = sup_g { E_Q[g] − log E_P[e^g] }`; the sup over `g` may be restricted to bounded continuous functions; proof in Dupuis & Ellis (1997), cited as Lemma 1.4.3 by secondary sources — [arXiv 2605.27253](https://arxiv.org/pdf/2605.27253); [arXiv 2310.17891](https://arxiv.org/pdf/2310.17891); [Shalizi lecture 31 (CMU 36-754)](https://www.stat.cmu.edu/~cshalizi/754/2006/notes/lecture-31.pdf)
- Rescaled "Gibbs tilt" form used in stochastic programming: `(1/γ) log E_P[e^{γX}] = sup_{Q ≪ P} { E_Q[X] − (1/γ) KL(Q‖P) }`, with the supremum attained by the Gibbs tilt `dQ*/dP ∝ e^{γX}` — [InFlow: entropic SDDP, arXiv 2609.11339](https://arxiv.org/pdf/2609.11339)
- Statistical-mechanics form (principle of minimal free energy): the Lagrangian `F = D − T H` (expected distortion minus temperature × Shannon entropy) is the Helmholtz free energy; its minimum over the association probabilities determines the equilibrium distribution, which is the Gibbs (canonical) distribution; the Lagrange multiplier is the temperature — [Rose 1998, Proc. IEEE, UCSB PDF](https://www.ece.ucsb.edu/publications/rose/pubs/pub42-Proc11-98.pdf)
- KL-divergence-constrained DRO (ambiguity in the objective) reduces to a one-layer convex minimization via the same duality; "KL-constrained DRO is often no harder than the underlying stochastic program" — [Hu & Hong, Optimization Online](https://optimization-online.org/?p=12225); Kocuk writes the same problem as a dual exponential-cone program solvable by conic solvers, with an uncapacitated facility location application — [Kocuk, IJOCTA 2021](https://www.accscience.com/journal/IJOCTA/11/2/10.11121/ijocta.01.2021.001001)
- In the SP literature the standard way to handle decision-dependent probabilities exactly is to make the probabilities variables and linearize products of binaries with continuous variables: Laumanns, Prestwich & Kawas show that neighbouring probability measures (in the space of influencing binary decisions) are linearly related via Bayes' rule, giving an exact MIP ("distribution shaping") — [HU Berlin edoc](https://edoc.hu-berlin.de/18452/9095); Basciftci, Ahmed & Shen linearize bilinear/trilinear terms of binary location variables times dual variables with McCormick envelopes (exact for the binary case, their Prop. 2) to get an exact MILP — [arXiv 1912.05577](https://arxiv.org/pdf/1912.05577)
- Hellemo, Tomasgard & Barton taxonomy: two-stage recourse models where prior probabilities are distorted by an affine transformation or blended as a convex combination of distributions, or where distribution parameters are first-stage decisions; distributions written exactly or via rational approximation; solved with BARON — [Springer, CMS 2018](https://link.springer.com/article/10.1007/s10287-018-0330-0); [MIT DSpace](https://dspace.mit.edu/handle/1721.1/117361)
- Bazotte, Carvalho & Vidal note that, unlike the exogenous case, convexity of the continuous relaxation of `Q(x, ξ)` does *not* yield convexity of the deterministic equivalent when probabilities depend on `x`, and that endogenous-probability DEPs are generally nonlinear nonconvex (NP-hard, citing Dupačová 2006) — [arXiv 2402.15486](https://arxiv.org/abs/2402.15486v3)

### Inferences
- **(EXACT) What DV gives directly.** Applying DV with `P = q`, `g_s = −β E_s(y) + t Q_s(y)`:
  `G_t(y) := (1/t) log Σ_s q_s exp(−β E_s(y) + t Q_s(y)) = max_{r ∈ Δ} { Σ_s r_s Q_s(y) − (β/t) Σ_s r_s E_s(y) − (1/t) KL(r‖q) }`.
  With `E_s` affine and `Q_s` convex, the inner function is convex in `y` (affine + convex, for fixed `r`) and concave in `r`, so `min_y G_t(y)` is a convex–concave saddle problem and `G_t` is convex in `y` (log-sum-exp composed with convex functions). But `G_t ≠ Σ_s p_s(y) Q_s(y)`; it equals `(1/t) log Z(y) + ρ_t^{p(y)}(Q(y))` where `ρ_t` is the entropic risk measure under `p(y)`.
- **(EXACT) The expectation is a derivative, not a sup.** `Σ_s p_s(y) Q_s(y) = ∂/∂t [ log Σ_s q_s exp(−β E_s(y) + t Q_s(y)) ]_{t=0}`. A derivative of a convex function is not convex in general, and there is no DV-type sup representation of an expectation by itself (the sup is attained only after exponentiation). This is why "y enters linearly in a min–max" cannot be obtained for the plain expectation.
- **(BOUND, DC) Sandwich by difference-of-convex functions.** Since `t ↦ log Z_t(y)` is convex, for any `t > 0`:
  `(1/t)[log Z(y) − log Z_{−t}(y)] ≤ Σ_s p_s(y) Q_s(y) ≤ (1/t)[log Z_t(y) − log Z(y)]`, where `Z_t(y) = Σ_s q_s exp(−β E_s(y) + t Q_s(y))`.
  Under affine `E` and convex `Q`, each `log Z_t(y)` is convex in `y`, so both bounds are explicit DC functions of `y` (convex minus convex); the gap is `O(t · max Var(Q))`. DC structure admits the convex–concave procedure (see KQ6).
- **(EXACT, derived here) Single-level mixed-integer convex reformulation via the Gibbs optimality gap.** For any `r ∈ Δ`, `KL(r‖p(y)) = Σ_s r_s log(r_s/q_s) + β Σ_s r_s E_s(y) + log Z(y) ≥ 0`, with equality iff `r = p(y)`. With affine energy, introduce `w_{sj} = r_s y_j` (McCormick: `w_{sj} ≤ r_s`, `w_{sj} ≤ y_j`, `w_{sj} ≥ r_s − 1 + y_j`, `w_{sj} ≥ 0`, exact when `y_j ∈ {0,1}`), and scale the recourse LP by `r_s` (`x_s := r_s x`, valid because the LP value is positively homogeneous in its right-hand side; `r_s = 0` forces `x_s = 0`). Then the original problem is equivalent to
  ```
  min   c_1^T y + Σ_s c^T x_s
  s.t.  W x_s ≥ r_s h_s − T_s w_s,  x_s ≥ 0                      ∀s
        McCormick(w_{sj}, r_s, y_j)                              ∀s,j
        Σ_s r_s = 1,  r ≥ 0
        Σ_s r_s log(r_s/q_s) + β Σ_s (a_s^T w_s + b_s r_s) + log Σ_s q_s exp(−β(a_s^T y + b_s)) ≤ 0
        y ∈ {0,1}^n ∩ Y.
  ```
  The last constraint is jointly convex in `(r, w, y)` (relative entropy + linear + log-sum-exp of affine), so the **continuous relaxation is a convex (exponential-cone representable) program** and the problem is a mixed-integer convex program solvable by outer approximation / branch-and-bound. `y` enters only linearly (in `w` and in the log-sum-exp argument). Caveats: (i) the relaxation is not tight because McCormick relaxes `w = r ⊗ y`; a cheap strengthening valid at every binary point is `Σ_s w_{sj} = y_j` (disjunctive hull of `{y_j=0, w=0} ∪ {y_j=1, w=r}`); (ii) **Slater's condition fails** for the convex constraint (its minimum value is exactly 0), so interior-point solvers may behave poorly; relaxing to `≤ ε` restores strict feasibility but makes the model optimistic (lower bound) since `r` is minimized over; (iii) exactness needs affine `E_s`. If `E_s` is merely convex, `log Z(y)` is not convex; if `E_s` is concave, `Σ_s r_s E_s(y)` is not convex — affine is the boundary of the exact-MICP regime.
- **(EXACT) Affine probabilities (Hellemo "affine distortion") are already MILP-representable** with binary `y`: `Σ_s (q_s + d_s^T y) Q_s(y)` needs only `y_j Q_s(y) = min{c^T x : W x ≥ y_j(h_s − T_s y)}` plus `y_j y_k` linearization; the Gibbs model is the natural *normalized, nonnegative* alternative, and the formulation above is its analogue.
- **(EXACT) The exact gradient contains a covariance ("score") term**: `∇_y Σ_s p_s(y) Q_s(y) = E_{p(y)}[∇Q(y)] − β Cov_{p(y)}(Q(y), a)`. Any scheme that "fixes p and optimizes y" drops the covariance term; this is only justified (by the envelope theorem) when the objective is the free energy, in which `p` is the argmin (KQ6, KQ8).

### Gaps
- I did not find the KL-gap MICP formulation above stated anywhere in the SP literature; it is derived here and unverified computationally. Its relaxation strength relative to the Laumanns et al. "distribution shaping" MIP or to Hellemo et al.'s BARON approach is unknown.
- Dupuis & Ellis's exact lemma number (1.4.3 vs 1.4.2) could not be verified against the book; only secondary citations were available.

---

## KQ2. Free energy and the log-partition function: with θ affine in y, which quantities are convex in y, and can the expected recourse be bounded/represented through them?

### Takeaway
With natural parameter `θ = −β y` and sufficient statistics `a_s`, `log Z(y)` is **exactly convex** in `y` and its gradient is `−β E_{p(y)}[a]` (so mean statistics are *monotone* in `y`, not convex). The expected recourse is a directional derivative of the cumulant generating function, hence not convex; it is representable exactly through `log Z` only by the KL-gap constraint of KQ1 or by DC bounds, and log-partition functions of *tilted* families give convex risk-sensitive surrogates.

### Cited Findings
- Conjugate duality of the log-partition function: for `μ` in the relative interior of the marginal polytope `MARG(G;φ)`, the conjugate dual `Φ*(μ)` equals the negative entropy of `p(x; θ(μ))` with `θ(μ)` dually coupled to `μ`; for `μ ∉ cl MARG`, `Φ*(μ) = +∞`; and `Φ(θ) = max_{μ ∈ MARG} {⟨θ, μ⟩ − Φ*(μ)}`, with the optimum attained uniquely at the exact marginals — [Wainwright & Jordan, NeurIPS 2003](https://proceedings.neurips.cc/paper_files/paper/2003/file/3cf2559725a9fdfa602ec8c887440f32-Paper.pdf); full treatment in the 2008 monograph, Sections 3.4 (marginal polytopes) and 3.6 (conjugate duality) — [Wainwright & Jordan 2008 PDF (Princeton mirror)](https://www.cs.princeton.edu/courses/archive/fall11/cos597C/reading/WainwrightJordan2008.pdf)
- Chernoff-type bounds on exponential-family probabilities can be posed as variational problems on the marginal polytope ("Variational Chernoff bounds for graphical models") — [arXiv 1207.4172](https://arxiv.org/pdf/1207.4172)
- Entropic risk / KL-DRO duality: `(1/γ) log E_P e^{γX} = sup_Q {E_Q X − (1/γ) KL(Q‖P)}` — [arXiv 2609.11339](https://arxiv.org/pdf/2609.11339); KL-DRO objective reduces to a one-layer convex minimization — [Hu & Hong](https://optimization-online.org/?p=12225)

### Inferences
- **(EXACT)** For `A(θ) = log Σ_s q_s exp(θ^T a_s)`: `A` convex, `∇A(θ) = E_θ[a]`, `∇²A(θ) = Cov_θ(a) ⪰ 0`. With `θ = −βy + const`: `log Z(y)` convex in `y`; `m(y) := E_{p(y)}[a] = −(1/β) ∇_y log Z(y)` is the gradient of a convex function, hence a monotone map (`(y − y')^T (m(y) − m(y')) ≤ 0`) but **not convex or concave componentwise**. Any objective of the form `λ^T m(y)` (recourse affine in the sufficient statistic and independent of `y`) is a directional derivative of a convex function — not convex.
- **(EXACT)** Quantities convex in `y`: `log Z(y)`; every tilted log-partition `log Σ_s q_s exp(−β E_s(y) + t Q_s(y))` for `t ≥ 0` (convex `Q_s`); the free energy `−(1/β) log Z(y)` is **concave**; `(1/t) log Z_t(y) − (1/t) log Z(y)` is DC.
- **(EXACT)** Representation through conjugate duality: `log Z(y) = sup_{μ ∈ conv{a_s}} { −β y^T μ − A*_q(μ) }` where `A*_q(μ) = min{KL(r‖q) : Σ_s r_s a_s = μ}` (relative-entropy I-projection value). So `−log Z(y) = inf_μ { β y^T μ + A*_q(μ) }` is a min over a fixed set of a jointly convex function — this is the free energy, and the "min–min" form the brief asks about holds for the *free energy*, not for the expected recourse.
- **(BOUND)** Jensen/Hölder: for `t > 0`, `Σ_s p_s(y) Q_s(y) ≤ (1/t) log Σ_s p_s(y) e^{t Q_s(y)} = (1/t)[log Z_t(y) − log Z(y)]`; the right side is DC in `y`. Equality iff `Q` is `p(y)`-a.s. constant.

### Gaps
- No source was found that bounds a decision-dependent expectation `E_{p(y)}[Q(y)]` by log-partition functions in `y`; the DC sandwich is derived here.

---

## KQ3. Mean-field and Bethe approximations: tractable surrogates exact for tree-structured or product-form p(y)?

### Takeaway
Mean-field gives a **lower bound on log Z** and Bethe a stationary-point approximation (exact on trees); tree-reweighted (TRW) gives a convex **upper bound**. These bound the *normalizer*, which in the finite-scenario setting is already cheap; they do not approximate `E_{p(y)}[Q(y)]` because the recourse LP couples all demand coordinates, so the only structures where they help are product-form demand with *separable* recourse, or as a replacement for `log Z(y)` inside the KQ1 constraint (which breaks exactness in a known direction).

### Cited Findings
- Mean-field and Bethe/loopy-BP are variational approximations of the log-partition function in the Wainwright–Jordan framework; mean field restricts the optimization to tractable (product) subfamilies and yields lower bounds, while Bethe replaces the exact entropy by the Bethe entropy, exact on trees — [Wainwright & Jordan 2008 PDF](https://www.cs.princeton.edu/courses/archive/fall11/cos597C/reading/WainwrightJordan2008.pdf) (Sections 5 and 4 per the monograph structure); the dual function is +∞ outside the closure of the marginal polytope, so relaxations must stay inside an outer bound of it — [Wainwright & Jordan 2003](https://proceedings.neurips.cc/paper_files/paper/2003/file/3cf2559725a9fdfa602ec8c887440f32-Paper.pdf)
- Rose's EM interpretation: viewing the association probability as the expectation of a binary assignment variable, the two-step DA iteration (compute Gibbs probabilities, then minimize over codevectors) is the EM algorithm; standard "zero-temperature" methods are the hard-assignment limit — [Rose 1998](https://www.ece.ucsb.edu/publications/rose/pubs/pub42-Proc11-98.pdf)
- Yuille & Rangarajan: all EM algorithms, Legendre-minimization and variational-bounding algorithms, and many mean-field theory algorithms are instances of CCCP — [Neural Computation 2003 (mlanthology)](https://mlanthology.org/neco/2003/yuille2003neco-concaveconvex)

### Inferences
- **(EXACT, product form)** If scenarios are joint demand vectors `d = (d_i)_i` and `p(d|y) ∝ q(d) exp(−β Σ_i E_i(d_i, y))` with `q` product-form, then `Z(y) = Π_i Z_i(y)`, `p(d|y) = Π_i p_i(d_i|y)`, and the KL-gap constraint of KQ1 can be written with per-customer marginals `r_i` (polynomial size). But `E_r[Q(y)]` requires the *joint* `r`, so the exactness of product/tree entropy does not transfer to the objective unless `Q_s(y) = Σ_i Q_i(d_i, y)` (separable recourse, e.g., independent lost-sales penalties with no shared capacity).
- **(APPROXIMATION, direction known)** Replacing `log Z(y)` in the KQ1 constraint by a mean-field *lower* bound loosens the constraint (the feasible set grows, giving an optimistic relaxation); replacing it by a TRW *upper* bound tightens it (may cut off the true `p(y)`, giving infeasibility or a pessimistic restriction). Bethe has no sign guarantee.
- **(APPROXIMATION)** Mean-field over `y` itself (optimization-as-inference, KQ6) replaces `f(y) = Σ_s p_s(y) Q_s(y)` by its multilinear extension `E_{y∼Bern(μ)}[f(y)]`, which is where products `y_i y_j` become `μ_i μ_j`; this is nonconvex in `μ` in general and is the "mean-field annealing" view of the bilinear facility-location set, not a convexification.

### Gaps
- No paper was found that applies Bethe/TRW bounds to a decision-dependent expectation in a stochastic program; the assessment is structural reasoning only.

---

## KQ4. Fenchel duality of log Z and the marginal polytope: a bridge to the convex hull of the bilinear facility-location set?

### Takeaway
In the scenario model the relevant "marginal polytope" is `conv{a_s}` (scenario statistics), not `conv{0,1}^n`; the conjugate of `log Z` is the I-projection value function. The genuine bridge to the facility-location bilinear set is the lifted variable `w = r ⊗ y`: the Gibbs/I-projection constraint lives in `(y, r, w)` space, and tightening the McCormick hull of `{(y, r, r⊗y) : y binary, r ∈ Δ}` is exactly the bilinear-hull question the OR literature already studies.

### Cited Findings
- `Φ*(μ) = −H(p(x;θ(μ)))` on the relative interior of the marginal polytope, `+∞` outside its closure; the marginal polytope for finite state spaces is a finite intersection of halfspaces — [Wainwright & Jordan 2003](https://proceedings.neurips.cc/paper_files/paper/2003/file/3cf2559725a9fdfa602ec8c887440f32-Paper.pdf)
- Basciftci, Ahmed & Shen model demand moments as piecewise-linear functions of `y` (`μ_j(y) = min{μ̄_j(1 + Σ_{j'} λ^μ_{jj'} y_{j'}), μ^UB_j}`, `σ_j^2(y) = max{σ̄_j^2(1 − Σ_{j'} λ^σ_{jj'} y_{j'}), (σ^LB_j)^2}`), use a moment ambiguity set `U(y)`, dualize, then linearize bilinear/trilinear products with McCormick envelopes (exact for products of binaries with bounded continuous variables, their Prop. 2) to obtain an exact MILP, plus valid inequalities — [arXiv 1912.05577](https://arxiv.org/pdf/1912.05577); [EJOR 2021 record](https://ideas.repec.org/a/eee/ejores/v292y2021i2p548-561.html)
- Luo & Mehrotra: decision-dependent ambiguity sets based on moment bounds, Wasserstein, φ-divergence and K–S test; finite-support reformulations via duality, solved by global optimization; their φ-divergence set has a fixed nominal `P_0` and decision-dependent radius `η(x)`, reformulated in Thm 4 as a program with constraint `α φ'(N p_i) + β + h(x, ξ^i) + λ_i = 0` — [arXiv 1806.09215](https://arxiv.org/pdf/1806.09215)
- Noyan, Rudolf & Lejeune: DRO with balls (earth mover's distances, including total variation and Wasserstein) centred at a **decision-dependent** probability distribution; applications to machine scheduling and humanitarian logistics — [IJOC 2022](https://ideas.repec.org/a/inm/orijoc/v34y2022i2p729-751.html)
- Maximum-capture facility location under MNL: captured demand is concave in `y`, enabling exact (multicut) outer approximation; Ljubić & Moreno combine outer-approximation and submodular cuts in branch-and-cut — [ESSEC record, EJOR 2018](https://faculty.essec.edu/research/outer-approximation-and-submodular-cuts-for-maximum-capture-facility-location-problems-with-random-utilities); [arXiv 2110.08497](https://arxiv.org/pdf/2110.08497)

### Inferences
- **(EXACT)** For the scenario family, `M = conv{a_s : s}` and `A*_q(μ) = min{KL(r‖q) : Σ_s r_s a_s = μ}` (convex, = `+∞` outside `M`). The decision `y` appears as the natural parameter `θ = −βy`; the marginal polytope of *binary y* (`[0,1]^n` with pairwise/higher marginals) would arise only if `y` were itself the random variable (optimization-as-inference, KQ6).
- **(EXACT, derived here)** The convex hull of `{(y_j, r, w) : w = r y_j, y_j ∈ {0,1}, r ∈ Δ}` for a single binary coordinate is `{ w ≥ 0, r − w ≥ 0, Σ_s w_s = y_j, Σ_s r_s = 1 }` (Balas disjunction of the two faces). This is strictly tighter than componentwise McCormick and costs one equality per facility. The hull over all `n` coordinates simultaneously is the hull of `Δ × {0,1}^n` lifted by all products — the same object as the bilinear facility-location set, and I found no polynomial-size description.
- **(EXACT)** The MNL concavity result is a special case where probabilities are *linear-fractional* in `y` (`p_{ij}(y) = y_j e^{u_ij} / (e^{u_i0} + Σ_k y_k e^{u_ik})`) and the "recourse" is trivial (captured weight). With a nontrivial `Q_s(y)`, the product `p_s(y) Q_s(y)` loses concavity; with binary `y` one can still use Charnes–Cooper/Dinkelbach plus exact linearization of `y_j Q_s(y)`, but this is parametric MILP, not a convex program.

### Gaps
- Whether the KL-gap formulation's convex relaxation dominates Basciftci et al.'s McCormick-based MILP relaxation (different probability models, so not directly comparable) is unknown.

---

## KQ5. Large deviations / Varadhan's lemma: use in risk-sensitive or decision-dependent stochastic programs?

### Takeaway
The Legendre-transform structure (`log E e^{nf} ≈ n sup(f − I)`) is exactly the DV duality; its operational use is the equivalence between **KL/φ-divergence DRO and entropic risk**, which is convex in `y` only when the centre distribution is fixed. With a decision-dependent Gibbs centre, the worst-case expectation is a DC function of `y` with explicit decomposition, and (derived here) is again a mixed-integer convex program after McCormick linearization for binary `y`.

### Cited Findings
- `KL(Q‖P) = sup_g {E_Q g − log E_P e^g}` with the sup over bounded measurable (or bounded continuous) functions; Varadhan–Donsker references — [arXiv 2605.27253](https://arxiv.org/pdf/2605.27253); [Shalizi lecture 31](https://www.stat.cmu.edu/~cshalizi/754/2006/notes/lecture-31.pdf)
- Hu & Hong: KL-constrained DRO ⇒ one-layer convex program (objective ambiguity); ambiguous chance constraints become the original chance constraint at an adjusted confidence level — [Optimization Online](https://optimization-online.org/?p=12225); exponential-cone formulation and facility-location application — [Kocuk 2021](https://www.accscience.com/journal/IJOCTA/11/2/10.11121/ijocta.01.2021.001001)
- Luo & Mehrotra Theorem 4: `D^3RO` with `P^φ(x) = {P : D_φ(P‖P_0) ≤ η(x)}` reformulated as a (semi-infinite or KKT-based) program; the inner problem is `max_p Σ_i p_i h(x, ξ^i)` s.t. `(1/N) Σ_i φ(N p_i) ≤ η(x)` — [arXiv 1806.09215](https://arxiv.org/pdf/1806.09215)
- Decision-dependent centre: Noyan, Rudolf & Lejeune use balls centred on a decision-dependent distribution (EMD-type) — [IJOC 2022](https://ideas.repec.org/a/inm/orijoc/v34y2022i2p729-751.html); Doan applies DRO under endogenous uncertainty to retrofitting planning — [Warwick WRAP](https://wrap.warwick.ac.uk/154930/1/WRAP-distributionally-robust-optimization-endogenous-uncertainty-application-retrofitting-planning-Doan-2021.pdf)
- Entropic SDDP uses the Gibbs-tilt form of the DV identity to build a risk-averse (entropic) multistage program — [arXiv 2609.11339](https://arxiv.org/pdf/2609.11339)

### Inferences
- **(EXACT)** KL ball around the Gibbs centre: `sup_{KL(r‖p(y)) ≤ ε} Σ_s r_s Q_s(y) = inf_{λ>0} { λ ε + λ log Σ_s q_s exp(−β E_s(y) + Q_s(y)/λ) − λ log Z(y) }`. For fixed `λ` and affine `E`, the middle term is convex in `y`, the last is concave in `y`: explicit DC decomposition; CCCP applies. As `ε → 0` this recovers `Σ_s p_s(y) Q_s(y)` (so the DC structure is inherited by the expectation only in the limit, where the `λ` scaling degenerates).
- **(EXACT, derived here)** For binary `y` and bounded `λ`, the whole worst-case functional is mixed-integer convex: write `−λ log Z(y) = min_{r} λ[KL(r‖q) + β E_r[E(y)]]`, substitute `ρ = λ r` (relative entropy `Σ_s ρ_s log(ρ_s/(λ q_s))` is jointly convex), linearize `ρ_s y_j` and `λ y_j` by McCormick, and note `λ log Σ_s q_s exp((−βλ E_s(y) + Q_s(y))/λ)` is the perspective of log-sum-exp applied to arguments that are affine in the lifted variables plus the convex `Q_s(y)` (epigraph `t_s ≥ c^T x_s`, `W x_s ≥ h_s − T_s y`), hence jointly convex. Needs a valid box on `λ` (a Slater/ε > 0 condition supplies one).
- **(EXACT)** Varadhan's lemma itself (`lim (1/n) log E e^{nF(X_n)} = sup_x {F(x) − I(x)}`) adds nothing beyond DV for finite scenario sets; its role is to justify entropic/exponential-utility objectives as the large-`n` limit of rare-event-weighted costs. No paper applying LDP asymptotics to decision-dependent two-stage programs was found.

### Gaps
- No source found that studies KL/φ-divergence ambiguity with a Gibbs/log-linear decision-dependent centre; Luo–Mehrotra make only the radius decision-dependent, Noyan et al. use EMD balls.

---

## KQ6. Variational inference view (Toussaint 2009; Levine 2018): ELBO as surrogate; does EM-style alternation give convergent coordinate descent with a convex y-subproblem?

### Takeaway
Control-as-inference yields a Boltzmann/soft policy and an ELBO equal to reward plus entropy, but Levine's key warning is that letting the variational distribution alter the *dynamics* makes the agent optimistic (risk-seeking) — the exact analogue of letting scenario probabilities `r` be free in a min–min. The "fix q, optimize y / fix y, optimize q" alternation is a monotone coordinate descent **only on the free-energy functional** (Rose's EM interpretation), where the `y`-subproblem is a convex fixed-probability stochastic program; applied to the expectation it is a heuristic fixed-point iteration that ignores the score term. The convergent, convex-subproblem version for the expectation is CCCP on the DC bound.

### Cited Findings
- Levine 2018: optimality variable `p(O_t = 1 | s_t, a_t) = exp(r(s_t, a_t))` (Sec. 2.2, Eq. 3; requires non-positive rewards, Sec. 2.5); variational `q(τ) = q(s_1) Π_t q(s_{t+1}|s_t,a_t) q(a_t|s_t)` with `q(s_1) = p(s_1)` and `q(s_{t+1}|s_t,a_t) = p(s_{t+1}|s_t,a_t)` fixed, so only the policy varies and dynamics terms cancel (Sec. 3.2, Eq. 18); bound `log p(O_{1:T}) ≥ E_q[ Σ_t r(s_t,a_t) − log q(a_t|s_t) ]` (Sec. 3.2, Eq. 19); optimal policy `π(a|s) = exp(Q(s,a) − V(s))`, `V(s) = log ∫ exp(Q(s,a)) da` (Sec. 3.1); exact inference with stochastic dynamics gives the backup `Q = r + log E[exp(V(s'))]` which "produces Q-functions that are optimistic" (Sec. 2.3, Eq. 7), and the posterior-dynamics objective yields "an optimistic policy that assumes a degree of control over the dynamics" (Sec. 2.4); fixing the dynamics resolves this (Sec. 3) — [ar5iv 1805.00909](https://ar5iv.labs.arxiv.org/html/1805.00909); abstract: max-entropy RL "is equivalent to exact probabilistic inference in the case of deterministic dynamics" — [arXiv 1805.00909](https://arxiv.org/abs/1805.00909)
- Toussaint 2009: a probabilistic model whose most likely trajectory coincides with the optimal one; reproduces LQG classically and extends to nonlinear systems via EP-like approximate inference — [ICML 2009 PDF](https://icml.cc/imls/conferences/2009/papers/271.pdf); [mlanthology record](https://mlanthology.org/icml/2009/toussaint2009icml-robot)
- Rose: the DA two-step iteration (Gibbs probabilities ↔ codevector update) is the EM algorithm; in the mass-constrained variant the free energy is minimized under the obvious constraints — [Rose 1998](https://www.ece.ucsb.edu/publications/rose/pubs/pub42-Proc11-98.pdf)
- CCCP: for an energy `E = E_vex + E_cave`, the iteration `∇E_vex(x^{t+1}) = −∇E_cave(x^t)` monotonically decreases `E`; all EM algorithms and variational-bounding algorithms can be re-expressed as CCCP — [NeurIPS 2001](https://papers.nips.cc/paper/2125-the-concave-convex-procedure-cccp); [Neural Computation 2003](https://mlanthology.org/neco/2003/yuille2003neco-concaveconvex)

### Inferences
- **(EXACT)** Levine's optimism pathology maps directly: `min_y min_{r ∈ Δ} Σ_s r_s Q_s(y)` without a KL penalty/constraint picks the cheapest scenario. The KL-gap constraint of KQ1 (or the penalty `KL(r‖p(y))`) is the SP analogue of "fixing the dynamics".
- **(EXACT)** Coordinate descent `r^{k+1} = argmin_r F(y^k, r)`, `y^{k+1} = argmin_y F(y, r^{k+1})` on `F(y, r) = Σ_s r_s Q_s(y) + (1/β) KL(r‖q) + Σ_s r_s E_s(y)` is monotone and the `y`-step is a convex fixed-probability stochastic program (`Q` convex, `E` affine); the `r`-step is closed-form Gibbs `r_s ∝ q_s exp(−β E_s(y) − β Q_s(y))`. But `min_r F = −(1/β) log Σ_s q_s exp(−β(E_s(y) + Q_s(y)))`, a soft-min/free energy, **not** `Σ_s p_s(y) Q_s(y)`. Iterating "set `r = p(y^k)`; minimize `Σ_s r_s Q_s(y)` over `y`" is not coordinate descent on any single function (because `p(y^k)` does not minimize `Σ_s r_s Q_s(y^k)` over `r`) and has no monotonicity guarantee; it drops the covariance term in the exact gradient (KQ1).
- **(BOUND, convergent)** CCCP on the DC upper bound `U_t(y) = (1/t)[log Z_t(y) − log Z(y)]`: at `y^k`, replace the concave part `−(1/t) log Z(y)` by its tangent `(β/t) m(y^k)^T y + const` with `m(y^k) = E_{p(y^k)}[a]`, and solve the convex (exponential-cone) subproblem `min_y c_1^T y + (1/t) log Σ_s q_s exp(−β E_s(y) + t Q_s(y)) + (β/t) m(y^k)^T y`. This is "fix q" entering through the mean sufficient statistic, is monotone in `U_t`, and converges to a stationary point of `U_t` (Yuille–Rangarajan); `U_t → Σ_s p_s(y) Q_s(y)` as `t → 0` but the subproblem's `1/t` scaling becomes ill-conditioned. Integrality of `y` must be handled by relaxing to `[0,1]^n` inside CCCP or by using CCCP as a primal heuristic.

### Gaps
- No work was found applying control-as-inference or ELBO surrogates to two-stage stochastic programs with decision-dependent scenario probabilities; the mapping is by analogy.

---

## KQ7. Entropic regularization / Sinkhorn: decision-dependent distribution as the entropic-OT (I-)projection of a reference onto a y-dependent constraint set?

### Takeaway
Entropic OT is a KL projection of a Gibbs kernel onto fixed marginal constraints; its decision-dependent analogue `p(y) = argmin{KL(r‖q) : A r = m(y)}` has a convex *value function* `A*_q(m(y))`, so the single-level constraint `KL(r‖q) ≤ A*_q(m(y))` is DC — unless one introduces the dual (natural) parameter `θ` as a variable, in which case the constraint `KL(r‖q) + log Z(θ) − θ^T m(y) ≤ 0` is convex after McCormick linearization of `θ ⊗ y` (binary `y`, bounded `θ`). The expectation becomes linear in the "plan" `r` only via the same `w = r ⊗ y` lifting; the Sinkhorn fixed-point form `diag(u) K diag(v)` is a nonconvex equality and is the wrong representation to put into an optimization model.

### Cited Findings
- Entropic OT solves `min ⟨P, M⟩ + reg · KL(P | K)` s.t. `P 1 = a`, `P^T 1 = b` with `K = exp(−M/reg)`; Sinkhorn iteratively projects onto the marginal constraints — [POT documentation](https://tessl.io/registry/tessl/pypi-pot/0.9.0/files/docs/entropic-transport.md)
- Peyré's numerical tour: entropic transport is the KL projection of the Gibbs kernel `K_{ij} = e^{−C_{ij}/ε}` onto couplings with prescribed marginals; the optimal plan has the form `diag(u) K diag(v)` and Sinkhorn alternates `u ← a/(Kv)`, `v ← b/(K^T u)` — [Numerical Tours notebook](https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/optimaltransp_5_entropic.ipynb)
- Sinkhorn is block optimization of a smooth convex function whose dual is a Bregman projection, with exponentially fast convergence — [arXiv 1812.11104](https://arxiv.org/pdf/1812.11104)
- Rose: earlier structured-DA work appealed to the principle of minimum cross-entropy (minimum divergence), a generalization of maximum entropy; alternate minimization of cross-entropy and leaf-layer cost "is not ensured to converge" — [Rose 1998](https://www.ece.ucsb.edu/publications/rose/pubs/pub42-Proc11-98.pdf)
- Luo (2021) derives decision-dependent demand from a "maximum attraction principle" (ranking of opened service centres, `q_ij ≤ y_j`), decoupled from the ambiguity set, giving MISOCP reformulations — [arXiv 2011.12514](https://arxiv.org/pdf/2011.12514)

### Inferences
- **(EXACT)** I-projection onto a linear family with `y`-dependent moments: `p(y) = argmin{KL(r‖q) : Σ_s r_s a_s = m(y)}` has the exponential-family form `p_s(y) ∝ q_s exp(θ(y)^T a_s)` with `θ(y)` defined implicitly (not affine in `y` even if `m(y)` is). Value function `A*_q(m) = sup_θ {θ^T m − log Σ_s q_s e^{θ^T a_s}}` is convex in `m`.
- **(EXACT, derived here)** Single-level reformulation: `min_{y,r,θ} Σ_s r_s Q_s(y)` s.t. `Σ_s r_s a_s = m(y)`, `KL(r‖q) + log Σ_s q_s e^{θ^T a_s} − θ^T m(y) ≤ 0`, `r ∈ Δ`, `y` binary. The last constraint forces `r = p(y)` and `θ = θ(y)`. With `m(y)` affine, `θ^T m(y)` is bilinear; for binary `y` and a box on `θ` (valid when `m(y)` stays in the relative interior of `conv{a_s}` uniformly), McCormick is exact and the constraint is convex (exponential cone). The objective is handled by the same `r_s`-scaling of the recourse LP as in KQ1. Contrast with KQ1: there `y` sits in the lower-level *objective* and the value function is concave (so the gap constraint is convex outright); here `y` sits in the lower-level *constraints*, the value function is convex, and convexity is recovered only by dualizing.
- **(EXACT)** Transport-plan form: if `p(y)` is a coupling `π(s, j)` between scenarios and facilities with marginals depending on `y` (e.g., `Σ_s π(s,j) = κ_j y_j`), then `Σ_{s,j} π(s,j) Q_s(y)` is linear in `π` only for fixed `y`; the product `π × Q(y)` needs the `w = π ⊗ y` lifting, so "expectation linear in a transport plan" holds after lifting, not before.

### Gaps
- No paper was found that models decision-dependent scenario probabilities as an entropic-OT/Sinkhorn projection inside a stochastic program.

---

## KQ8. Physics-inspired work on facility location / stochastic programming (deterministic annealing, free-energy objectives): do they handle decision-dependent probabilities?

### Takeaway
Deterministic annealing (Rose; Salapaka and co-workers) *does* optimize an expectation under a decision-dependent Gibbs distribution — the association probabilities `p(j|i; y) ∝ exp(−d(x_i, y_j)/T)` depend on the facility decisions — but the energy *is* the cost, so the free energy `F = D − T H` is a legitimate smoothed surrogate recovered exactly at `T → 0`, and the alternation is EM on `F`. None of this work treats an exogenous scenario distribution tilted by decisions with a separate recourse cost; for that setting the transferable ideas are the KL-gap constraint (KQ1), DC/CCCP (KQ6), and annealing in `β` or `t` as a homotopy.

### Cited Findings
- Rose 1998 (Proc. IEEE 86:2210–2239): DA minimizes the Lagrangian `F = D − T H` (Helmholtz free energy) where `D` is expected distortion and `H` Shannon entropy; the minimizing association probabilities are Gibbs distributions parameterized by `T`; annealing tracks the free-energy minimum while lowering `T`; the system undergoes phase transitions at computable critical temperatures; at `T = 0` assignments are hard and standard methods (GLA/k-means) are the zero-temperature limit; the two-step iteration is EM; mass-constrained clustering needs only as many codevectors as distinct clusters — [Rose 1998 PDF](https://www.ece.ucsb.edu/publications/rose/pubs/pub42-Proc11-98.pdf); background on Rose's thesis (association probabilities at a given average distortion are Gibbs, Lagrange multiplier = inverse temperature) — [Caltech thesis record](https://thesis.caltech.edu/2858/)
- Salapaka group: DA for clustering/facility location with capacity and size constraints, minimizing a free-energy function over association probabilities with inverse temperature as annealing parameter — [arXiv 1606.06427](https://arxiv.org/pdf/1606.06427); entropy-based handling of inequality constraints in facility location — [arXiv 2002.03505](https://ar5iv.labs.arxiv.org/html/2002.03505); capacitated facility location with utilization bounds in a maximum-entropy framework using control barrier functions (2025) — [arXiv 2504.01378](https://arxiv.org/html/2504.01378v1); free energy as a coverage measure with resource dynamics keeping its time derivative non-positive — [CDC-ECC 2005](https://skoge.folk.ntnu.no/prost/proceedings/cdc-ecc05/pdffiles/papers/2762.pdf); [ACC 2008](https://folk.ntnu.no/skoge/prost/proceedings/acc08/data/papers/1192.pdf); DA for multiple TSP — [arXiv 1604.04169](https://arxiv.org/pdf/1604.04169)
- Decision-dependent SP solution approaches that are *not* physics-based but exact: distribution shaping (Bayes'-rule linear relations ⇒ MIP, Istanbul highway instance solved to optimality in seconds) — [HU Berlin edoc](https://edoc.hu-berlin.de/18452/9095); random-variable transformation mapping endogenous to exogenous uncertainty, with Bernoulli (inversion), exponential (inversion), normal (standardization), binomial (convolution of Bernoullis) handled, yielding MILP or mixed-integer convex programs under stated criteria, solved by SAA — [arXiv 2402.15486](https://arxiv.org/abs/2402.15486v3); Hellemo et al. (affine distortion / convex combination / parameters as decisions; BARON) — [Springer CMS 2018](https://link.springer.com/article/10.1007/s10287-018-0330-0); scenario aggregation with BDDs for endogenous uncertainty — [arXiv 1701.04055](https://arxiv.org/pdf/1701.04055)
- MNL-based (Gibbs-type) facility location: concave objective enables exact outer approximation; submodular cuts — [ESSEC record](https://faculty.essec.edu/research/outer-approximation-and-submodular-cuts-for-maximum-capture-facility-location-problems-with-random-utilities); maximum-capture under random utilities with capacities/overflow penalties (2025) — [UC Chile repository](https://repositorio.uc.cl/dspace/items/e05bf1ac-f43a-4356-879d-b4ed93f8a1f9)

### Inferences
- **(EXACT)** In DA, `F(y) = −T Σ_i log Σ_j exp(−d(x_i, y_j)/T) = min_{p} { D(y, p) − T H(p) }` and `D(y) = Σ_{i,j} p(j|i; y) d(x_i, y_j) = F(y) + T H(p(y))`, so `F ≤ D ≤ F + T Σ_i log K` (`K` = number of facilities). The envelope theorem gives `∇_y F = E_{p(y)}[∇_y d]`, i.e. the DA/EM `y`-step correctly ignores the dependence of `p` on `y` **because `p(y)` is the argmin of `F`**. This is the precise sense in which "free energy objectives" make decision-dependent Gibbs probabilities harmless — and it requires the energy to equal the cost.
- **(EXACT)** For the user's problem the analogous free energy is `−(1/β) log Σ_s q_s exp(−β(E_s(y) + Q_s(y)))` (KQ6), a risk-seeking soft-min of `E + Q`, not the risk-neutral expectation under `p(y)`. Hence DA is not directly applicable; what transfers is the homotopy idea (anneal `t` in the DC bound `U_t`, or `β`) and the EM/CCCP structure.
- **(EXACT)** `F(y)` in DA is not convex in `y` (log-sum-exp of concave functions), which is why annealing/phase transitions are needed; convexity in the user's setting comes only from affine energy plus the lifted `(r, w)` variables (KQ1), not from the free-energy form itself.

### Gaps
- No physics-inspired (DA, mean-field annealing, free-energy) paper was found that treats a two-stage stochastic program whose scenario probabilities are an exogenous distribution tilted by first-stage decisions; the closest are MNL-capture models (trivial recourse) and DA facility location (energy = cost).
- I did not open Salapaka et al.'s full texts to confirm whether any variant separates a demand-scenario distribution from the assignment energy; the abstracts indicate they do not.

---

## Summary table for the operations researcher

| Identity | Status | Requires | What it yields for `min_y Σ_s p_s(y) Q_s(y)` | Convex in `y`? |
|---|---|---|---|---|
| DV: `log Σ q_s e^{f_s} = max_r {E_r f − KL(r‖q)}` | EXACT | bounded `f` | representation of tilted log-partition / entropic risk, not of the expectation | the tilted functional `G_t` is convex (affine `E`, convex `Q`); the expectation is not |
| Gibbs gap: `KL(r‖p(y)) = KL(r‖q) + β E_r[E(y)] + log Z(y)` | EXACT | affine `E_s`, binary `y` | single-level mixed-integer convex (exp-cone) program with lifted `w = r⊗y`, scaled recourse `x_s = r_s x` | relaxation convex; Slater fails; relaxation not tight |
| `log Z(θ)` convex, `∇ log Z = E[a]` | EXACT | exponential family | `log Z(y)` convex; `E_{p(y)}[a]` monotone; expectation = ∂_t of convex function | no |
| DC sandwich `(1/t)[log Z_t − log Z]` | BOUND | affine `E`, convex `Q` | explicit DC bounds with `O(t)` gap; CCCP (convex subproblems) | DC only |
| Mean-field / Bethe / TRW | BOUND / APPROX | graphical scenario structure | bounds on `log Z(y)` only; help objective only for separable recourse | n/a |
| Conjugate `A*(μ) = −H` on marginal polytope | EXACT | finite scenarios | `conv{a_s}` is the relevant polytope; bridge to FL bilinear hull is via `w = r⊗y` | hull strengthening `Σ_s w_{sj} = y_j` |
| KL-DRO ⇔ entropic risk (Varadhan/DV) | EXACT | fixed centre ⇒ convex; Gibbs centre ⇒ DC / MICP after lifting | decision-dependent KL-DRO with explicit DC decomposition | DC; MICP for binary `y` |
| Control-as-inference ELBO / EM alternation | EXACT for free energy | energy = cost | monotone coordinate descent on free energy; convex `y`-step; optimistic if `r` free | no (only subproblems) |
| Entropic OT / I-projection with `y`-dependent moments | EXACT | affine `m(y)`, bounded dual `θ`, binary `y` | MICP after dualizing and McCormick on `θ⊗y`; DC without dualization | convex after lifting |
| Deterministic annealing (Rose; Salapaka) | EXACT for `F`, `F ≤ D ≤ F + T log K` | energy = cost | homotopy in `T`; not applicable when scenario energy ≠ recourse cost | no |
