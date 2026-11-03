# notes/

Research notes behind the paper. Reading order:

1. `how-the-method-was-found.md` — the path from the question to the solver: survey, verification, the per-factor
   lifted cuts, the linking rows, the scaling walls, the scenario decomposition, the C++ solver, the referee rounds,
   and the claims that were withdrawn.
2. `reframing-bilinearity-report.md` — the cross-field survey (2026-10-28): every reformulation device for a
   decision-dependent expectation, ranked, with sources. Its first table is the map of the field.
3. `reframing-bilinearity/` — the seven field surveys the report was written from, and the verification note:
   - `or_baseline.md` — exact methods in stochastic programming (Pantuso & Hewitt, Basciftci et al., Laumanns, Noyan).
   - `global_opt_submodularity.md` — RLT, perspective hulls, multilinear polytopes; the origin of the lifted cuts.
   - `variational_physics.md` — Gibbs/KL, Donsker-Varadhan, DC sandwiches (Candidate 1).
   - `economics_bilevel.md` — Charnes-Cooper, choice models, bilevel collapse (Candidate 3).
   - `mdp_rl_ideas.md` — occupation measures, performative convexity, Hammer-Rubin (Candidate 4).
   - `change_of_measure.md` — likelihood ratios, push-in couplings (the convex mean-shift alternative).
   - `experimental_design_information.md` — value of information, posterior lifts.
   - `verification_of_derived_reformulations.md` — proofs and toy-instance checks of the four candidates; the
     per-factor cuts (2b) and the level-1 RLT LP (2c) are identified here.
4. `literature-instance-sizes.md` — largest instances solved to optimality in the comparable literature, read from
   the PDFs (2026-10-31), used in the related-work section.

The experiment record is `experiments/RESULTS.md`; the paper is `paper/main.tex`.
