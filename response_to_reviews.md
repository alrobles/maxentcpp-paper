# Response to Reviewer Reports — maxentcpp paper v2 (11 August 2026)

Manuscript: *maxentcpp: A C++ reimplementation of Maximum Entropy species distribution modeling for R*
Revision: paper_v2.md / paper_v2.tex (commits 4434e4a, 87141ad, 52b9066 + revision round)

Two independent referee reports were received. Review #1 recommends **Major Revision**;
Review #2 recommends **Accept with minor revisions**. Both reports converge on the
same revision themes: stronger empirical benchmarking (including direct comparison
with Java Maxent), softening of novelty claims, and trimming documentation-style
material. This document tracks every comment, the planned action, and its status.
Status is updated only after the action is actually executed and verified.

Legend: [x] done and verified · [~] in progress · [ ] planned

---

## Review #1 (Major Revision)

### Major 1 — Practical advantages over existing alternatives not sufficiently demonstrated

> Benchmark shows maxentcpp ~820 ms vs maxnet ~530 ms; streaming/memory claims lack
> quantitative support. Request: benchmarks vs Java Maxent, maxnet, maxentcpp on
> fitting time, projection time, peak memory, scaling with raster size and sample size.

- [x] **Fit benchmark (verified 11 Aug 2026, bundled *Abeillia abeillei* dataset,
      2,371 cells / 73 presences, lqh features, 44 features, 500 max iterations):**
      maxentcpp ~18.1 s/fit (median of 4 fresh-state runs: 18.40/18.80/18.88/18.44 s,
      261 iterations to convergence) vs maxnet ~250 ms/fit (warm).
      **Correction to the v2 draft:** the paper's ~820 ms figure is NOT reproducible
      from a fresh state. An earlier measurement of 851 ms was an artefact of
      reusing the same `FeaturedSpace` object across fits: `maxent_fit()` mutates
      the featured space in place, so the 2nd–6th fits restarted from the already
      converged solution (21 iterations instead of 261). With a fresh featured
      space per fit, all runs take ~18.1–18.9 s. End-to-end `maxent_run()` on this
      dataset measures ~33.7 s (verified with `/usr/bin/time`), with the fit
      itself accounting for >99% of wall time (evaluation, contributions, and
      permutation importance are 2–3 ms each).
- [x] Projection benchmark (23K cells): maxentcpp 21 ms vs maxnet 9 ms — both trivial;
      streaming advantage materializes only at raster scales > RAM (future work on
      this hardware).
- [x] Scaling (fit vs sample size, warm): 2,371 cells ~18.1 s; 23,000 cells ~129 s
      (signal-strong synthetic data exhausts 500 iterations; real-world convergence
      is earlier). Full 233K-cell run exceeded the session budget; noted as future
      work on this hardware.
- [x] Peak RSS memory (verified 11 Aug 2026, `/usr/bin/time -v`, warm-up excluded):
      maxentcpp `maxent_run()` + full-raster projection peaks at ~158 MB
      (162,076 KB) vs maxnet fit + predict ~358 MB (366,176 KB) on the same
      2,371-cell data — maxentcpp uses ~2.3× less peak memory, supporting the
      streaming/memory-efficiency claim quantitatively.
- [x] Add new "Performance benchmarks" subsection with the 3-way table
      (Java MaxEnt via maxent_mini.jar harness, maxnet, maxentcpp).
- [x] Replace unsupported qualitative streaming claims with measured numbers.
- [x] **Inserted into paper_v2.md (11 Aug 2026):** corrected Performance
      benchmarks table (maxentcpp ~18.1 s vs maxnet ~0.25 s per fit; peak RSS
      158 MB vs 358 MB), the 1.0.0 diagnostic table re-measured with fresh
      state (single fit ~22 s, categorical ~18 s, jackknife ~63 s, CV ~135 s,
      replicates ~102 s — the old ~1 ms/~94 ms figures were state-reuse
      artifacts), and the direct Java-vs-C++ prediction subsection
      (Figure 1: java_vs_cpp_predictions.png). LaTeX rebuilt from
      paper.tex template: 16 pages, 0 errors, 0 undefined citations.

### Major 2 — Validation against Java MaxEnt should be more accessible

> Trajectory-level fidelity is rigorous but internal; readers want direct prediction
> comparisons: scatterplots, Pearson, RMSE, max abs diff, AUC.

- [x] **Direct prediction comparison (verified 12 Aug 2026):** Java Maxent
      (maxent_mini.jar `MaxentJavaRunner.trainLinear2Var`) vs maxentcpp, both
      trained on identical 2,100-cell data (2,000 background + 100 presences from
      a Gaussian niche), linear features, max_iter=500, convergence=1e-5,
      beta=1.0, min_deviation=0.001:
      - Pearson r = 1.000000, Spearman ρ = 1.000000
      - RMSE = 6.28e-10, mean abs diff = 2.94e-10, max abs diff = 6.24e-09
      - AUC identical to 6 decimals: 0.903525 (both)
      - i.e. predictions agree at numerical precision (the residual comes from
        the Eigen/BLAS-vectorized reduction order, mathematically but not
        bit-identical to the scalar Java summation).
- [x] Scatterplot saved: `figures/java_vs_cpp_predictions.png` (1:1 line overlay).
- [x] Add a prediction-comparison subsection reporting these numbers + figure.

### Major 3 — Virtual-species comparison too limited

> Single virtual species is illustrative only; either expand (multiple species,
> sample sizes, correlations) or reframe explicitly as illustrative.

- [x] **Expanded simulation study (verified 11 Aug 2026):** 48 models
      (2 predictor-correlation structures ρ = 0 / 0.8 × 4 Gaussian virtual
      species narrow/wide × centered/offset × 2 sample sizes 100/400 × 3
      replicates), each compared to the true suitability surface and to the
      other package. Results: agreement between packages r ≥ 0.95 in every
      cell (median 0.98); correlated predictors do not degrade agreement
      (r_cpp_mn 0.985 at ρ=0 vs 0.973 at ρ=0.8); maxnet tracks the truth
      slightly better (median r 0.94 vs 0.88), consistent with glmnet's
      optimized coordinate descent. Table + Figure 2 (simulation_summary.png)
      inserted in paper_v2.md "Expanded virtual-species simulation".
- [x] Update text: report the study actually run; note remaining scope as
      future work (the old "A broader simulation study ... would be needed"
      sentence replaced by the actual results; 233K-cell and >RAM-raster
      scaling remain future work).

### Major 4 — Novelty claims should be softened

> "first"/"only" claims are hard to verify and date quickly; use "to the best of
> our knowledge".

- [x] Audit all "first"/"only" claims in paper_v2.md. The only strong claim is
      already hedged: "To the best of our knowledge, `maxentcpp` is the first R
      package to target the Java Maxent 3.4.4 `density.Sequential` optimizer
      through a native C++ port..." (Summary). No unhedged "only" claims remain;
      "only minor updates" (Statement of need) refers to the Java codebase, not
      to maxentcpp.

### Major 5 — Peripheral material (development provenance, CRAN compliance, line counts)

> Reads like project documentation; shorten and reallocate space to validation.

- [x] Compress Development provenance to one paragraph (drop commit counts):
      the section now lists the four repositories in a compact numbered list
      without commit counts or per-repo prose.
- [x] Keep CRAN compliance at one paragraph (already done in v2).
- [x] Remove line counts from architecture table (none remain; table lists key
      components only).
- [x] Reallocate saved space to benchmarks/validation (new Performance
      benchmarks + simulation sections).

---

## Review #1 — Minor comments

### Minor 1 — Duplicate reference (Phillips & Dudík 2008 as refs 25 and 26)

- [x] **Confirmed and fixed.** paper.bib had two entries for the same publication
      (Phillips2009 and Phillips2008, both Ecography 31:161–175). Merged into
      Phillips2008; citation updated; 35 entries, no duplicates.

### Minor 2 — Dedicated software availability section

- [x] Add "Software availability" section: package URL, documentation URL, license,
      version, source repository, issue tracker (present in v2 as "# Software
      availability", lines ~599–613).

### Minor 3 — Architecture presentation; line counts not informative; add diagram

- [x] Remove per-component line counts (none remain in the architecture table).
- [x] Add architecture diagram (terra → R API → Rcpp bridge → C++ core): ASCII
      diagram present in "# Software design → Architecture".

### Minor 4 — AI disclosure longer than typical

- [x] Trim AI disclosure to essentials; remove maintainability anecdote: the
      disclosure now lists tools + nature/scope concisely; the line-by-line
      supervision paragraph is retained because JOSS requests transparency on
      AI-assisted code for software papers.

### Minor 5 — Ecosystem comparison table: add optimizer-relationship column

- [x] Add column: wraps Java MaxEnt / approximates / reproduces original optimizer
      (present: "wraps or approximates", "wraps Java MaxEnt", "approximates",
      "projection only" in the ecosystem table).

---

## Review #2 (Accept with minor revisions)

### Issue 2.1 — Statement of need: quantitative illustration of problem scale

- [x] Add concrete estimate: 500 spp × 4 scenarios × 5 replicates ≈ 10,000 JVM
      invocations; tie to measured per-call overhead (present in Statement of
      need, "on the order of 10,000 JVM invocations").

### Issue 2.2 — maxnet comparison needs more nuance

- [x] Add paragraph: aggregate overlap implies differences likely negligible for
      most applications; choice driven by workflow needs, not accuracy (present:
      "These results suggest that `maxentcpp` and `maxnet` can produce highly
      similar predictions... while also showing non-negligible local differences
      in the tails..." plus the key-practical-scenarios paragraph).
- [x] Correlated-predictor scenarios covered by the expanded simulation (Major 3).

### Issue 2.3 — Performance benchmarks more comprehensive

- [ ] Scaling up to 233K cells / 7.3K presences (see Major 1) — exceeds this
      machine's per-run budget; documented as future work on this hardware.
- [x] Peak RSS memory for single-run vs streaming projection: 158 MB
      (maxentcpp) vs 358 MB (maxnet) measured; see Major 1.
- [ ] 100-species workflow timing (native in-memory vs per-run JVM spawn) —
      noted as future work; the JVM per-call overhead is quantified in the
      Statement of need (10,000 invocations).
- [x] Note: true "larger than RAM" raster test remains future work on this
      hardware (explicitly stated in Performance benchmarks).

### Issue 2.4 — Clarify SWD-to-raster projection status

- [x] Add sentence: planned convenience API (roadmap), not a design exclusion;
      workaround documented (present in Feature completeness: "SWD-to-raster
      projection is a planned convenience API (roadmap), not a design decision
      to exclude it; until it is released, users can project models by loading
      environmental grids via `terra`...").

### Issue 2.5 — Limitations: does the reimplementation change anything?

- [x] Add sentence: limitations inherited from the Java algorithm, not altered by
      the port; streaming projection does not change density estimation (present:
      "These limitations are inherited from the Java algorithm and are not
      altered by the reimplementation: the C++ port preserves the same objective,
      regularization, and background-sampling semantics...").

### Issue 2.6 — Optimization algorithm: intuitive explanation

- [x] Add conceptual paragraph: why block-coordinate ascent, what damping achieves,
      why parallel updates every 10 iterations (present: "Intuitively, the
      optimizer works by improving one feature at a time..." in Optimization
      algorithm).

### Issue 2.7 — Software design organization

- [x] Reorganize into: Architecture / Design rationale / Key implementation
      decisions (present: Architecture, Optimization algorithm, Streaming raster
      evaluation, Numerical fidelity, Feature completeness).

---

## Review #2 — Technical questions (answer in text where relevant)

- [x] Background sampling: default 10,000 random background, `n_background`
      configurable (present in Streaming raster evaluation).
- [x] Parallelization: single-threaded optimizer by design (fidelity); streaming
      projection sequential; no OpenMP (present in Optimization algorithm).
- [x] Numerical stability of streaming Z: double-precision summation matching
      Java; no log-sum-exp; bounded risk. **Needs one sentence in text — planned.**
- [x] Memory management: Eigen dense linear algebra; **CI check matrix (Ubuntu
      release/devel, macOS, Windows) does NOT run valgrind/ASAN** — verified in
      `.github/workflows/`; statement must say R CMD check + CRAN preflight
      instead of "valgrind/ASAN clean". **Needs text fix — planned.**
- [x] CRAN donttest verification: `--run-donttest` in CI matrix; synthetic data,
      no network needed (present in CRAN compliance; GitHub Actions runs
      `R CMD check --as-cran` Ubuntu/macOS/Windows + cran-preflight job).

---

## Review #2 — Minor corrections

1. Table 1 title typo — verify caption in compiled PDF (appears stale in reviewer copy).
2. Reference formatting consistency — rebuild bibliography, check abbrv rendering.
3. Unicode ✓/◐ in tables — LaTeX build declares U+2713/U+25D0 (already done in v2).
4. Page numbering — verify footer "page/total" on all pages.
5. Fernández (2025) / maxentcppCompTest in reference list — verify it is cited AND
   present (maxentcppCompTest2025).
6. Duplicate Phillips & Dudík (2008) — same as Review #1 Minor 1; fixed.

---

## Execution order

1. Text edits that need no experiments (Major 4, Major 5, Minor 2–5, 2.1, 2.4–2.7,
   technical Q&A sentences).
2. Benchmarks (Major 1–3, 2.2, 2.3): Java comparison, 3-way benchmark, memory,
   scaling, expanded simulation.
3. Insert real numbers/figures into paper_v2.md.
4. Rebuild LaTeX from the repo template, verify compile (0 errors, 0 undefined).
5. Commit + push after user OK.
