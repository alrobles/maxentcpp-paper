---
title: 'maxentcpp: A C++ reimplementation of Maximum Entropy species distribution modeling for R'
tags:
  - R
  - C++
  - species distribution modeling
  - maximum entropy
  - ecology
  - biogeography
  - Rcpp
authors:
  - name: Ángel Luis Robles Fernández
    orcid: 0000-0002-4741-8012
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Ecology and Evolutionary Biology, University of Kansas, Lawrence, KS, United States
    index: 1
date: 11 May 2026
bibliography: paper.bib
---

# Summary

`maxentcpp` is an R package offering a native C++17
reimplementation of the Maximum Entropy (Maxent) algorithm for species
distribution modeling (SDM). Maxent estimates the geographic distribution
of a species from presence-only occurrence records and environmental
covariates by identifying the probability distribution of maximum entropy
subject to constraints derived from training data
[@Phillips2006; @Phillips2004]. Since its introduction, Maxent and its
foundational methodological papers have been cited extensively across
ecology, conservation biology, and biogeography [@Elith2011].

The fitted model is a *Gibbs distribution* over geographic space:

$$
P(\mathbf{x}) = \frac{1}{Z}\exp\!\left(\sum_{j=1}^{J} \lambda_j f_j(\mathbf{x})\right),
$$

where $f_j$ are feature transformations of the environmental variables,
$\lambda_j$ are learned weights, and
$Z = \sum_{\mathbf{x}}\exp\!\bigl(\sum_j \lambda_j f_j(\mathbf{x})\bigr)$
is the normalizing constant (partition function). The optimizer finds the
$\lambda_j$ that maximize the regularized log-likelihood via sequential
coordinate ascent with $\ell_1$ (lasso) penalties
[@Phillips2006; @Elith2011].

The original Maxent software is a Java desktop application
[@Phillips2017]. `maxentcpp` translates the core continuous-predictor
Maxent 3.4.4 fitting and projection algorithm into C++17 with R bindings
via Rcpp [@Eddelbuettel2011] and RcppEigen [@Bates2013], eliminating the
Java Runtime Environment (JRE) dependency entirely. The package supports all five feature types (linear,
quadratic, product, hinge, and threshold), output transformations in raw,
logistic, and complementary log-log (cloglog) scales, and a streaming
raster evaluation engine that processes environmental layers block-by-block
through the `terra` package [@Hijmans2024], enabling continental-scale
projections on commodity hardware without loading entire raster stacks
into memory. Tools for diagnosing model performance include AUC
evaluation, permutation importance, response curves, and Multivariate
Environmental Similarity Surfaces (MESS) for novelty detection.

# Statement of need

Maxent is the *de facto* standard for presence-only species distribution
modeling [@Merow2013]. However, the original Java implementation presents
a number of practical barriers for present-day ecological research
workflows:

1. **Java dependency.** Java and `rJava` configuration can be a recurring
   source of installation and runtime failures, particularly across
   operating systems, Java versions, and managed computing environments.
   While `conda` environments and Docker containers can manage Java
   versioning, they add deployment complexity and are not standard
   practice in many ecology research groups.

2. **File-based I/O.** Data must be serialized to disk as SWD or ASCII grid
   files before the Java process can read them, and results must be parsed
   back from output files. There is no in-memory bridge between R data
   structures and the Java optimizer.

3. **Scalability.** The Java implementation is single-threaded and
   GUI-centric. In multi-species, multi-scenario workflows (e.g., 500
   species $\times$ 4 climatic scenarios), file I/O overhead accumulates:
   each species run requires writing SWD files, invoking the JVM, and
   parsing output files. A native in-memory API eliminates this per-call
   overhead entirely.

4. **Maintenance.** The Java Maxent codebase has received only minor
   updates since version 3.4.4, with limited active development
   [@Phillips2017].

`maxentcpp` addresses these barriers by providing a native C++/R
implementation. At present, the development version can be installed from
GitHub using `remotes::install_github("alrobles/maxentcpp")`; following
CRAN acceptance, the package will be installable with
`install.packages("maxentcpp")`. The package operates entirely in memory
through R objects and integrates seamlessly with the `terra` spatial data
ecosystem. The intended users are ecologists, conservation biologists,
and biogeographers who employ Maxent in reproducible, script-based
workflows — from individual species analyses to large-scale biodiversity
assessments.

# Legacy software modernization in ecology

Computational ecology relies on a small number of foundational software
tools, many of which were written over a decade ago and have not been
modernized for current computing environments. Reimplementing legacy
scientific software in modern languages is uncommon in ecology but has
important precedents in adjacent environmental sciences.

`maxentcpp` follows the same modernization logic as projects such as
Circuitscape.jl [@Hall2021]: remove brittle runtime dependencies, improve
integration with present-day workflows, and retain scientific
comparability with established software.

Closer to species distribution modeling, `maxnet` [@Phillips2024maxnet]
represents the original MaxEnt author's own reimplementation using R and
`glmnet` --- preserving the same statistical model but employing a
different optimization algorithm (coordinate descent on the full feature
matrix rather than sequential coordinate ascent).

`maxentcpp` advances this modernization trend by pursuing a stricter
compatibility objective than most wrapper packages: preserving the Java
Maxent `density.Sequential` optimizer trajectory rather than solely
reproducing final predictions. To the best of current knowledge,
`maxentcpp` is the first R package to target the Java Maxent 3.4.4
`density.Sequential` optimizer through a native C++ port and to publish
per-iteration trajectory tests against the Java implementation via the
companion package `maxentcppCompTest` [@maxentcppCompTest2025].

# State of the field

Several R packages provide access to Maxent-family models.
**dismo** [@Hijmans2023] is the most widely used R interface to Java
Maxent; it wraps `maxent.jar` via `rJava`, inheriting the JDK dependency
and file I/O overhead, and is currently in maintenance mode following the
transition from `raster` to `terra`.
**maxnet** [@Phillips2024maxnet] is a pure-R implementation that uses
`glmnet` for elastic-net regularization; while it avoids the Java
dependency, it employs a different optimization algorithm (elastic net via
coordinate descent on the full feature matrix) rather than the sequential
coordinate-ascent optimizer of the original Maxent, which can produce
numerically different results, particularly for threshold and hinge
features.
**ENMeval** [@Kass2021] is a model tuning and evaluation framework that
wraps either `dismo::maxent()` or `maxnet::maxnet()`; it does not
implement the Maxent algorithm itself.
**kuenm** [@Cobos2019] is a calibration and evaluation toolkit that
relies on `dismo` or direct Java Maxent calls.

## Empirical comparison: maxentcpp vs maxnet

To quantify the practical differences between `maxentcpp` and `maxnet`,
we compared both packages on a virtual species [@Leroy2016] with known
true habitat suitability, constructed from the environmental layers
bundled with `maxentcpp` (Annual Mean Temperature and Annual
Precipitation over Central America, 2,371 non-NA cells). The virtual
species has a Gaussian niche centered at 20 °C and 1,500 mm, and 100
presence records were sampled proportionally to the true suitability
surface
(see the companion vignette
[Virtual Species Comparison](https://alrobles.github.io/maxentcpp/articles/virtual_species_comparison.html)
for full reproducible code).

Both packages were fitted with linear, quadratic, and hinge features.
In this illustrative example, the two packages produced highly similar
rank and spatial-overlap patterns:

| Metric | Value |
|--------|------:|
| Schoener's $D$ | 0.93 |
| Spearman $\rho$ (rank correlation) | 0.95 |
| Pearson $r$ | 0.97 |
| Mean $|\Delta\text{cloglog}|$ | 0.053 |
| Max $|\Delta\text{cloglog}|$ | 0.43 |

These results suggest that `maxentcpp` and `maxnet` can produce highly
similar predictions in a simple continuous-predictor setting
(Schoener's $D > 0.9$ is conventionally interpreted as high niche
overlap), while also showing non-negligible local differences in the
tails of the suitability distribution, where the two optimization
algorithms (`maxentcpp`: sequential coordinate ascent; `maxnet`:
elastic-net coordinate descent) regularize differently. A broader
simulation study across multiple species, sample sizes, and predictor
correlation structures would be needed to establish general equivalence.

The key practical scenarios where a user must choose `maxentcpp` over
`maxnet` are: (1) when exact reproduction of Java Maxent results is
required (e.g., replicating published analyses), (2) when lambda file
interoperability with Java Maxent is needed, and (3) when streaming
raster projection onto grids larger than available RAM is required.

## Performance benchmarks

Training time for the bundled *Abeillia abeillei* dataset (73 presences,
2,371 background, linear + quadratic + hinge features):

| Package | Median time per fit |
|---------|--------------------:|
| `maxentcpp` | ~820 ms |
| `maxnet` | ~530 ms |

`maxnet` is faster for small datasets because `glmnet`'s coordinate
descent is highly optimized for dense feature matrices.
`maxentcpp`'s streaming evaluation avoids materializing the full
prediction matrix, which is expected to reduce peak memory during
projection; future benchmarks should quantify this advantage on rasters
larger than available RAM.

`maxentcpp` was built as a new package rather than contributing to
existing projects for three reasons. First, a faithful port of the
original Java optimizer (sequential coordinate ascent using
$\ell_1$-penalized Newton steps) requires a C++ engine that cannot be
expressed through `glmnet`'s API; contributing this to `maxnet` would
change its core algorithm. Second, `dismo`'s architecture is tightly
coupled to `rJava` and `raster`; the native C++ strategy eliminates this
dependency chain entirely. Third, the header-based C++ library design of
`maxentcpp` enables reuse in other packages or standalone applications
beyond R --- something not achievable through modifications to existing
R-only packages.

# Software design

## Architecture

`maxentcpp` is organized in three layers (C++ core, Rcpp bridge, R
interface), with the C++ core further split into algorithmic and
infrastructure sublayers:

| Layer | Key components | Lines |
|-------|----------------|------:|
| C++ core (algorithmic) | `Sequential`, `FeaturedSpace`, five feature types | ~4,500 |
| C++ core (I/O, diagnostics) | `BackgroundProvider`, grid I/O, CSV, MESS, response curves | ~2,400 |
| Rcpp bridge | `rcpp_*.cpp` external-pointer bindings [@Eddelbuettel2013] | ~3,300 |
| R interface | `maxent_run()`, feature generators, projection, evaluation, diagnostics | ~2,500 |

The C++ core uses Eigen [@Guennebaud2010] for dense linear algebra.
Of the ~6,900 C++ lines, approximately 4,500 implement the core
algorithmic logic (optimizer, features, density) while the remainder
handles I/O and diagnostics. The high-level
`maxent_run()` function provides a one-call workflow mirroring the Java
GUI experience, while lower-level functions
(`maxent_generate_features()`, `maxent_featured_space()`,
`maxent_fit()`, `maxent_project_cloglog()`) offer full control over each
modeling step.

## Optimization algorithm

The `Sequential` optimizer is a direct port of `density.Sequential` in
Java Maxent 3.4.4, as represented by the public `mrmaxent/Maxent`
repository and AMNH release notes (where 3.4.4 is described as a minor
bug-fix release). It minimizes the $\ell_1$-regularized negative
log-likelihood:

$$
\mathcal{L}(\boldsymbol{\lambda}) = -\frac{1}{m}\sum_{i=1}^{m}
\sum_{j=1}^{J} \lambda_j f_j(\mathbf{x}_i) + \log Z(\boldsymbol{\lambda})
+ \sum_{j=1}^{J} \beta_j |\lambda_j|,
$$

where $m$ is the number of presence records, $Z(\boldsymbol{\lambda})
= \sum_{k=1}^{N} \exp\!\bigl(\sum_j \lambda_j f_j(\mathbf{x}_k)\bigr)$
is the partition function summed over $N$ background points, and
$\beta_j = \hat\sigma_j / \sqrt{m}$ is the per-feature regularization
parameter with $\hat\sigma_j$ being the standard deviation of feature
$j$ over the presence points (floored at 0.001).

At each iteration, the optimizer selects the feature $j^*$ with the
most negative $\Delta\mathcal{L}$ bound (the `deltaLossBound` function),
then computes a Newton step:

$$
\alpha_j = -\frac{\partial \mathcal{L}/\partial \lambda_j}
{\mathbf{u}_j^\top \mathbf{H} \mathbf{u}_j},
$$

where $\mathbf{u}_j$ is the $j$-th coordinate direction and
$\mathbf{u}_j^\top \mathbf{H} \mathbf{u}_j = \mathrm{Var}_q[f_j]$
is the variance of feature $j$ under the current distribution $q$.
The derivative includes the $\ell_1$ subgradient:
$\partial \mathcal{L}/\partial \lambda_j = \mathbb{E}_q[f_j] -
\mathbb{E}_{\hat p}[f_j] + \beta_j \,\mathrm{sign}(\lambda_j)$.
The step is damped during early iterations
($\alpha \leftarrow \alpha/50$ for iterations $<10$,
$\alpha/10$ for $<20$, $\alpha/3$ for $<50$) and falls back to
a line search (`searchAlpha`) if the Newton step does not decrease loss.
Every 10 iterations, a parallel update applies coordinate steps to all
features simultaneously.

Convergence is tested every 20 iterations: training stops when the
relative loss improvement falls below $10^{-5}$ or 500 iterations are
reached. All constants in this section are taken directly from Java
Maxent 3.4.4 source files, and each statement is covered by a
source-mapping test in `maxentcppCompTest`.

## Streaming raster evaluation

`maxentcpp` reads `terra::SpatRaster` objects block-by-block through a
`BackgroundProvider` abstraction. Each tile is scored by the C++ engine
and discarded before the next tile is loaded, allowing projection onto
raster stacks larger than available RAM. The partition function $Z$ is
accumulated across tiles by summing unnormalized densities
$\exp(\sum_j \lambda_j f_j(\mathbf{x}_k))$ for each cell $k$ in the
tile, then normalizing once after all tiles are processed. This produces results equivalent to single-pass evaluation because the
mathematical sum is commutative and each cell is scored independently;
floating-point summation order may introduce differences at the
least-significant bit level.

## Numerical fidelity

Each C++ class is a direct translation of its Java counterpart,
preserving the same control flow, regularization logic, and numerical
constants:

| `maxentcpp` (C++) | Java Maxent original |
|-------------------|----------------------|
| `LinearFeature` | `density.LinearFeature` |
| `QuadraticFeature` | `density.SquareFeature` |
| `ProductFeature` | `density.ProductFeature` |
| `HingeFeature` | `density.HingeFeature` |
| `ThresholdFeature` | `density.ThresholdFeature` |
| `FeaturedSpace` | `density.FeaturedSpace` |
| `Sequential::run()` | `density.Sequential.run()` |

This design choice prioritizes backward compatibility: users migrating
from Java Maxent should obtain numerically equivalent results for
implemented features under matched inputs, defaults, and deterministic
settings. The fidelity contract is
enforced by a dedicated companion package, `maxentcppCompTest`
[@maxentcppCompTest2025], which runs the C++ optimizer and the original
Java `density.Sequential` on shared test fixtures and compares
per-iteration trajectories of loss, entropy, and $\lambda$ vectors.

On controlled test fixtures (identical input data, same floating-point
representation), agreement reaches $< 10^{-14}$ relative error for
symmetric fixtures and $< 10^{-6}$ on $\|\Delta\lambda\|_\infty$
for asymmetric fixtures at every iteration checkpoint.
**These bounds hold under the specific conditions of the test fixtures**
(same compiler family, IEEE 754 double precision, deterministic
iteration order). Floating-point ordering differences introduced by
alternative BLAS backends or aggressive compiler optimizations
(e.g., `-ffast-math`) could widen the gap, though the sequential
nature of the coordinate-ascent algorithm limits sensitivity to
parallelism-induced reordering. Cross-platform CI (Ubuntu, macOS,
Windows with GCC, Clang, and MSVC) confirms that the test fixtures
pass on all three compiler families.

**Lambda file compatibility.** Trained models are serialized in the same
`.lambdas` text format used by Java Maxent 3.4.4. Lambda files produced by
either implementation can be loaded by the other, ensuring
interoperability with existing workflows and archived models.

## Feature completeness

The following table summarizes the current scope of `maxentcpp` relative
to Java Maxent 3.4.4 and `maxnet`:

Legend: ✓ = implemented and tested; ◐ = partial or experimental;
— = not implemented.

| Feature | `maxentcpp` | Java Maxent 3.4.4 | `maxnet` 0.1.4 |
|---------|:-----------:|:-----------:|:--------:|
| Linear features | ✓ | ✓ | ✓ |
| Quadratic features | ✓ | ✓ | ✓ |
| Product features | ✓ | ✓ | ✓ |
| Hinge features | ✓ | ✓ | ✓ |
| Threshold features | ✓ | ✓ | ✓ |
| Raw/logistic/cloglog output | ✓ | ✓ | ✓ |
| AUC evaluation | ✓ | ✓ | via ENMeval |
| Permutation importance | ✓ | ✓ | via ENMeval |
| Response curves | ✓ | ✓ | ◐ |
| MESS maps | ✓ | ✓ | — |
| Clamping / fade-by-clamping | ✓ | ✓ | — |
| Lambda file I/O | ✓ | ✓ | — |
| Streaming raster projection | ✓ | — | — |
| Bias grids | ✓ | ✓ | via `offset` |
| Categorical variables | — | ✓ | ✓ |
| Replicate runs / cross-validation | — | ✓ | via ENMeval |
| Jackknife variable selection | — | ✓ | — |
| Missing data handling | — | ✓ | ✓ |
| SWD-to-raster projection | — | ✓ | — |

Categorical variables, replicate runs, jackknife, and SWD-to-raster
projection are not yet implemented. Users requiring these features should
continue using Java Maxent or `maxnet`. We plan to add categorical
variable support in a future release.

## Development provenance

The translation followed a phased approach across four publicly
accessible repositories:

1. `alrobles/Maxent` --- a fork of the original `mrmaxent/Maxent`
   [@Phillips2017] Java source code, where the architecture was analyzed
   and each Java class was mapped to a C++ target (193 commits).
2. `alrobles/maxentcpp-devel` --- the development repository where the C++
   port, Rcpp bindings, R interface, tests, and documentation were
   iteratively built and refined (34 commits).
3. `alrobles/maxentcppCompTest` --- the cross-language fidelity test suite
   containing Java oracle runners, golden trajectory files, and automated
   comparison scripts (40 commits).
4. `alrobles/maxentcpp` --- the release repository with the clean,
   CRAN-ready package (47 commits).

# Ecosystem comparison

The R ecosystem contains several packages that provide access to MaxEnt
models, each employing a distinct integration strategy:

| Package | MaxEnt integration | Dependency |
|---------|-------------------|------------|
| **dismo** [@Hijmans2023] | rJava bridge to `maxent.jar` | Java JDK, rJava |
| **kuenm** [@Cobos2019] | `system2("java -jar maxent.jar")` | Java JDK |
| **kuenm2** (Cobos et al.) | `glmnet` via forked `maxnet` code | glmnet |
| **wallace** [@Kass2018wallace] | Delegates to ENMeval $\to$ maxnet or dismo | ENMeval |
| **ENMTools** [@Warren2010] | `dismo::maxent()` directly | dismo, rJava |
| **biomod2** [@Thuiller2009] | Java (`system2`) or `maxnet` | Optional Java or maxnet |
| **rmaxent** [@Baumgartner2017] | Java-free projection of fitted Maxent models, `.lambdas` parsing, MESS | None (projection only) |
| **MIAmaxent** [@Vollering2019] | Maximum entropy via subset selection (not lasso); decoupled transformation/fitting/selection | glmnet |
| **SDMtune** [@Vignali2020] | Trains/tunes SDMs via `dismo::maxent()` or `maxnet`; hyperparameter tuning, variable selection | dismo or maxnet |
| **maxentcpp** | Native C++17 via Rcpp | None (self-contained) |

These packages can be grouped by their relationship to the MaxEnt
algorithm:

- **Black-box wrappers** (dismo, kuenm, ENMTools, biomod2 MAXENT mode):
  invoke the Java binary and read output files. The optimizer is
  inaccessible.
- **Approximate reimplementations** (maxnet, kuenm2, biomod2 MAXNET
  mode): use `glmnet` coordinate descent on the logistic regression
  formulation. The statistical model is equivalent but the optimization
  path differs, which can produce numerically different results for
  threshold and hinge features.
- **Faithful reimplementation** (maxentcpp): ports the original
  Sequential optimizer to C++ and proves per-iteration numerical
  equivalence.

The distinguishing contribution of `maxentcpp` is that it is the only
package offering a compiled native reimplementation of the actual MaxEnt
`density.Sequential` optimizer --- preserving both the statistical model
and the optimization algorithm while eliminating the Java dependency.

# Limitations and inherited assumptions

The Maxent 3.4.4 algorithm that `maxentcpp` faithfully reproduces has
well-documented limitations that users should be aware of:

- **Feature explosion.** The number of features (particularly hinge and
  threshold) grows with the number of environmental variables, which can
  cause identifiability issues in high-dimensional settings
  [@Elith2011; @Merow2013].
- **Sensitivity to background sampling.** Maxent's output is influenced
  by the choice and size of the background sample, which affects the
  estimated density and can bias predictions in undersampled regions
  [@Phillips2009].
- **$\ell_1$ regularization limitations.** The sequential coordinate
  ascent with $\ell_1$ penalties produces sparse models but does not
  guarantee selection of the "correct" features under high correlation
  among predictors [@Hastie2015].
- **No principled uncertainty quantification.** Unlike Bayesian
  approaches or bootstrap-based methods, the Maxent optimizer produces
  point estimates without confidence intervals.

A principled alternative would be to reformulate the Maxent objective
using modern convex optimization tooling (e.g., proximal gradient
methods, ADMM) or Bayesian inference. However, such reformulations
would produce different results than the widely used Java implementation,
breaking comparability with published analyses. The design of
`maxentcpp` explicitly aims to preserve this comparability.

# Research impact statement

`maxentcpp` currently provides a numerically faithful implementation of
the core continuous-predictor Maxent fitting and projection workflow for
the implemented feature classes and output transformations. It is not yet
a complete replacement for all Java Maxent workflows, because categorical
predictors, replicate runs, jackknife diagnostics, missing-data handling,
and SWD-to-raster projection are not implemented in the current release.
For those features, users should continue using Java Maxent or `maxnet`.

The package is distributed under the MIT license (compatible with
Eigen's MPL2 and RcppEigen's GPL-2 via the "or later" clause), with
full source code, a pkgdown documentation website at
<https://alrobles.github.io/maxentcpp/>, and four vignettes covering a
getting-started walkthrough, the mathematical foundations of Maxent
features, a comparative motivation document, and a virtual species
validation study. Continuous integration via GitHub Actions runs
`R CMD check` on Ubuntu, macOS, and Windows across R release and R-devel.
The test suite comprises over 20 test files (~4,800 lines) covering
feature generation, optimization convergence, spatial projection, Java
numerical equivalency, model diagnostics, and end-to-end workflows.

By providing a dependency-free, memory-efficient, and numerically faithful
Maxent implementation, `maxentcpp` may lower the barrier for large-scale
biodiversity assessments, climate-change projections, and conservation
prioritization studies that rely on presence-only species distribution
models, once the remaining Java Maxent workflow features --- especially
categorical predictors and replicate diagnostics --- are implemented or
clearly documented as out of scope.

# AI usage disclosure

Generative AI tools were used during the development of `maxentcpp` and
the preparation of this manuscript, as detailed below.

**Tools used.**

- **GitHub Copilot** (GPT-4-based, 2025 version): code scaffolding and
  boilerplate generation during the initial porting phases in
  `alrobles/Maxent` and `alrobles/maxentcpp-devel`.
- **Devin** (Cognition AI, 2025--2026): incremental feature
  implementation, test writing, documentation generation, CRAN compliance
  fixes, CI configuration, and manuscript drafting.
- **Claude** (Anthropic, 2025--2026): complex refactoring, cross-cutting
  documentation, and code review.

**Nature and scope of assistance.** AI tools contributed to: C++ code
generation and refactoring from Java source, Rcpp binding scaffolding,
R wrapper functions, `testthat` test scaffolding and expansion, `roxygen2`
documentation, vignette drafting, pkgdown site configuration, CI workflow
setup, CRAN compliance review, and manuscript drafting. The core
algorithmic C++ classes (`Sequential`, `FeaturedSpace`, feature types)
were translated from Java with AI assistance but under direct human
supervision: each class was mapped line-by-line from the corresponding
Java source, reviewed for correctness against the Java reference, and
validated through the `maxentcppCompTest` trajectory comparison framework.
Approximately 60% of the C++ algorithmic lines were initially drafted by
AI tools and subsequently reviewed and corrected by the human author;
the remaining 40% were written directly by the author, particularly
the performance-critical optimizer inner loops and numerical edge cases.

**Maintainability.** The human author (Á. L. Robles Fernández) can
independently explain and modify every algorithmically significant class.
As evidence: the `Sequential::run()` optimizer, `good_alpha()`,
`delta_loss_bound()`, `newton_step_feature()`, and `do_parallel_update()`
methods are documented with line-by-line references to the corresponding
Java source (e.g., "mirrors Sequential.java:294..342"), and the author
has independently debugged numerical discrepancies during the
fidelity testing process.

**Confirmation of review.** All AI-generated code, tests, documentation,
and manuscript text were reviewed, edited, and validated by the human
author, who made all core design decisions including the algorithm
translation strategy, API design, numerical fidelity requirements, and
the phased development architecture. The full development history is
publicly available in the four repositories listed above.

# Acknowledgements

The author thanks Steven J. Phillips, Robert P. Anderson, and Robert
E. Schapire for creating the original Maxent software and making the
source code publicly available under an MIT license. The `Rcpp`,
`RcppEigen`, and `terra` package maintainers are gratefully acknowledged
for providing the infrastructure that makes `maxentcpp` possible.

# References
