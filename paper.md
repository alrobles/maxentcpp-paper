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

`maxentcpp` is an R package that provides a high-performance C++17
reimplementation of the Maximum Entropy (Maxent) algorithm for species
distribution modeling (SDM). Maxent estimates the geographic distribution
of a species from presence-only occurrence records and environmental
covariates by finding the probability distribution of maximum entropy
subject to constraints derived from the training data
[@Phillips:2006; @Phillips:2004]. Since its introduction, Maxent has
become one of the most widely used methods in ecology and conservation
biology, with over 13,000 citations [@Elith:2011].

The original Maxent software is a Java desktop application
[@Phillips:2017]. `maxentcpp` translates the complete Maxent 3.4.4
algorithm into C++17 with R bindings via Rcpp [@Eddelbuettel:2011] and
RcppEigen [@Bates:2013], eliminating the Java Runtime Environment (JRE)
dependency entirely. The package supports all five feature types (linear,
quadratic, product, hinge, and threshold), the sequential coordinate-ascent
optimizer with L1 regularization, and output transformations in raw,
logistic, and complementary log-log (cloglog) scales. A streaming
raster evaluation engine processes environmental layers block-by-block
through the `terra` package [@Hijmans:2024], enabling continental-scale
projections on commodity hardware without loading entire raster stacks
into memory. Model diagnostics include AUC evaluation, permutation
importance, response curves, and Multivariate Environmental Similarity
Surfaces (MESS) for novelty detection.

# Statement of need

Maxent is the *de facto* standard for presence-only species distribution
modeling. However, the original Java implementation presents several
practical barriers for modern ecological research workflows:

1. **Java dependency**: Users must install and configure a compatible JDK.
   Version conflicts between Java 8, 11, 17, and 21 produce silent failures,
   and the `rJava` bridge package is one of the most common sources of
   installation errors in the R ecosystem, particularly on macOS and HPC
   clusters.

2. **File-based I/O**: Data must be serialized to disk as SWD or ASCII grid
   files before the Java process can read them, and results must be parsed
   back from output files. There is no in-memory bridge between R data
   structures and the Java optimizer.

3. **Scalability**: The Java implementation is single-threaded and
   GUI-centric. Large-scale studies involving thousands of species, multiple
   climate scenarios, or ensemble pipelines are bottlenecked by file I/O
   overhead and the inability to integrate Maxent as a library call within
   a larger R pipeline.

4. **Maintenance**: The Java Maxent codebase has remained essentially
   unchanged since version 3.4.4 (released November 2020), with only two
   active contributors.

`maxentcpp` addresses these barriers by providing a native C++/R
implementation that installs with a single `install.packages("maxentcpp")`
call, operates entirely in memory through R objects, and integrates
seamlessly with the `terra` spatial data ecosystem. The target audience is
ecologists, conservation biologists, and biogeographers who use Maxent in
reproducible, script-based workflows — from individual species analyses to
large-scale biodiversity assessments.

# State of the field

Several R packages provide access to Maxent-family models:

- **dismo** [@Hijmans:2023]: The most widely used R interface to Java Maxent.
  It wraps `maxent.jar` via `rJava`, inheriting the JDK dependency and
  file I/O overhead. The package is in maintenance mode following the
  transition from `raster` to `terra`.

- **maxnet** [@Phillips:2024]: A pure-R implementation of Maxent using
  `glmnet` for elastic-net regularization. While `maxnet` avoids the Java
  dependency, it uses a different optimization algorithm (elastic net via
  coordinate descent on the full feature matrix) rather than the sequential
  coordinate-ascent optimizer of the original Maxent. This can produce
  numerically different results, particularly for threshold and hinge
  features.

- **ENMeval** [@Kass:2021]: A model tuning and evaluation framework that
  wraps either `dismo::maxent()` or `maxnet::maxnet()`. It does not
  implement the Maxent algorithm itself.

- **kuenm** [@Cobos:2019]: A calibration and evaluation toolkit for
  ecological niche models that relies on `dismo` or direct Java Maxent
  calls.

`maxentcpp` was built as a new package rather than contributing to
existing projects for three reasons. First, a faithful port of the
original Java optimizer (sequential coordinate ascent with
$\ell_1$-penalized Newton steps) requires a C++ engine that cannot be
expressed through `glmnet`'s API; contributing this to `maxnet` would
change its core algorithm. Second, `dismo`'s architecture is tightly
coupled to `rJava` and `raster`; the native C++ approach eliminates this
dependency chain entirely. Third, the header-only C++ library design of
`maxentcpp` enables reuse in other packages or standalone applications
beyond R, something not achievable through modifications to existing
R-only packages.

# Software design

## Architecture

`maxentcpp` is organized in three layers:

1. **C++ core** (~10,000 lines across 18 header and 12 source files):
   A header-based library under `src/cpp/include/maxent/` implementing the
   Maxent algorithm independently of R. Key classes include `FeaturedSpace`
   (the Gibbs distribution manager and trainer), `Sequential` (the
   coordinate-ascent optimizer ported from Java `density.Sequential`), and
   the five feature types (`LinearFeature`, `QuadraticFeature`,
   `ProductFeature`, `HingeFeature`, `ThresholdFeature`). The C++ layer
   uses Eigen [@Guennebaud:2010] for vectorized linear algebra.

2. **Rcpp bridge** (~1,500 lines): Bindings in `src/rcpp_*.cpp` that
   expose C++ objects to R as external pointers, following the Rcpp
   modules pattern [@Eddelbuettel:2013].

3. **R interface** (~4,000 lines across 15 files): User-facing functions
   that handle input validation, spatial data conversion via `terra`, and
   output formatting. A high-level `maxent_run()` function provides a
   one-call workflow mirroring the Java GUI experience, while lower-level
   functions (`maxent_generate_features()`, `maxent_featured_space()`,
   `maxent_fit()`, `maxent_project_cloglog()`) offer full control over
   each modeling step.

## Key design decisions

**Streaming raster evaluation.** Rather than materializing the full
background matrix in memory, `maxentcpp` reads `terra::SpatRaster`
objects block-by-block through a `BackgroundProvider` abstraction. Each
tile is scored by the C++ engine and discarded before the next tile is
loaded. This allows projection onto raster stacks larger than available
RAM.

**Numerical fidelity over algorithmic novelty.** Every C++ class is a
direct port of the corresponding Java class, preserving the same
control flow, regularization logic, and numerical constants. This design
choice prioritizes backward compatibility: users migrating from Java
Maxent obtain identical results. The fidelity contract is enforced by a
dedicated companion package, `maxentcppCompTest`
[@maxentcppCompTest:2025], which runs the C++ optimizer and the original
Java `density.Sequential` on shared test fixtures and compares per-iteration
trajectories of loss, entropy, and lambda vectors. On standardized
fixtures, agreement reaches $< 10^{-14}$ relative error.

**Lambda file compatibility.** Trained models are serialized in the same
`.lambdas` text format used by Java Maxent 3.4.4. Lambda files produced by
either implementation can be loaded by the other, ensuring
interoperability with existing workflows and archived models.

## Development provenance

The translation followed a phased approach across four repositories:

1. `alrobles/Maxent` — a fork of the original `mrmaxent/Maxent`
   [@Phillips:2017] Java source code, where the architecture was analyzed
   and each Java class was mapped to a C++ target (193 commits).
2. `alrobles/maxentcpp-devel` — the development repository where the C++
   port, Rcpp bindings, R interface, tests, and documentation were
   iteratively built and refined (34 commits).
3. `alrobles/maxentcppCompTest` — the cross-language fidelity test suite
   containing Java oracle runners, golden trajectory files, and automated
   comparison scripts (40 commits).
4. `alrobles/maxentcpp` — the release repository with the clean,
   CRAN-ready package (47 commits).

# Research impact statement

`maxentcpp` is designed to serve as a drop-in replacement for Java Maxent
in R-based ecological workflows. The package eliminates the most common
installation barrier (JRE + `rJava` configuration), enables integration
of Maxent into fully scripted and reproducible pipelines, and supports
raster-scale projections that are infeasible with file-based I/O.

The package is distributed under the MIT license, with full source code,
documentation, and a pkgdown website at
<https://alrobles.github.io/maxentcpp/>. Continuous integration runs
`R CMD check` on Ubuntu, macOS, and Windows across R release and R-devel.
The test suite comprises over 20 test files covering feature generation,
optimization convergence, spatial projection, Java numerical equivalency,
model diagnostics, and end-to-end workflows.

By providing a dependency-free, memory-efficient, and numerically faithful
Maxent implementation, `maxentcpp` lowers the barrier for large-scale
biodiversity assessments, climate-change projections, and conservation
prioritization studies that rely on presence-only species distribution
models.

# AI usage disclosure

Generative AI tools were used during the development of `maxentcpp` and
the preparation of this manuscript, as follows:

**Software development.** GitHub Copilot (GPT-4-based, 2025 version) was
used for code scaffolding and boilerplate generation during the initial
porting phases in `alrobles/Maxent` and `alrobles/maxentcpp-devel`. Devin
(Cognition AI, 2025–2026) was used for incremental feature implementation,
test writing, documentation generation, CRAN compliance fixes, and CI
configuration across `alrobles/maxentcpp-devel`, `alrobles/maxentcpp`, and
`alrobles/maxentcppCompTest`. Claude (Anthropic, 2025–2026) was used for
complex refactoring, cross-cutting documentation, and code review.

**Paper authoring.** Devin (Cognition AI, 2026) assisted with drafting the
initial text of this manuscript based on the repository documentation,
commit history, and JOSS formatting requirements.

**Scope of assistance.** AI tools contributed to: C++ code generation and
refactoring from Java source, Rcpp binding scaffolding, R wrapper
functions, `testthat` test scaffolding and expansion, roxygen2
documentation, vignette drafting, pkgdown site configuration, CI workflow
setup, CRAN compliance review, and manuscript drafting.

**Human review.** All AI-generated code, tests, documentation, and
manuscript text were reviewed, edited, and validated by the human author
(Á. L. Robles Fernández), who made all core design decisions including
the algorithm translation strategy, API design, numerical fidelity
requirements, and the phased development architecture. The full
development history is publicly available in the four repositories listed
above.

# Acknowledgements

The author thanks Steven J. Phillips, Robert P. Anderson, and Robert
E. Schapire for creating the original Maxent software and making the
source code publicly available under an MIT license. The `Rcpp`,
`RcppEigen`, and `terra` package maintainers are gratefully acknowledged
for providing the infrastructure that makes `maxentcpp` possible.

# References
