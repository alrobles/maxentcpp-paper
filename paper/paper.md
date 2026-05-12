---
title: 'maxentcpp: A C++ reimplementation of Maximum Entropy species distribution modeling for R'
tags:
  - R
  - Rcpp
  - species distribution modeling
  - Maxent
  - ecological modeling
  - biogeography
  - conservation biology
  - raster
  - terra
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

`maxentcpp` is an R package providing a native C++17 implementation of the
core continuous-predictor Maximum Entropy (Maxent) species distribution
modeling workflow. It reproduces the fitting and projection algorithm of Java
Maxent 3.4.4 [@Phillips2017] entirely in C++ via Rcpp [@Eddelbuettel2011] and
RcppEigen [@Bates2013], eliminating the Java/JDK runtime dependency that has
been a persistent barrier in ecological modeling workflows. The package
integrates with `terra` [@Hijmans2024] for modern raster I/O and supports
linear, quadratic, product, hinge, and threshold feature classes with raw,
logistic, and complementary log-log (cloglog) output scales.

# Statement of need

Maxent is one of the most widely used algorithms for modeling species
geographic distributions from occurrence records and environmental covariates
[@Phillips2006; @Elith2011]. The foundational methodological papers have been
cited extensively across ecology, conservation biology, and biogeography.
Despite this widespread adoption, the original Java implementation introduces
practical barriers for present-day R-based workflows:

1. **JDK dependency**: Java Maxent requires a compatible Java Development Kit
   and the `rJava` bridge package, which frequently fails to configure
   correctly across operating systems, institutional computing clusters, and
   containerized environments.
2. **File-based I/O**: The Java implementation communicates through temporary
   files rather than in-memory data structures, adding overhead and
   complicating integration with modern spatial packages.
3. **Maintenance status**: The `dismo` package [@Hijmans2023], which provides
   the primary R interface to Java Maxent, is in maintenance mode following
   the transition from `raster` to `terra`.
4. **Reproducibility**: Java version conflicts between JDK 8, 11, 17, and 21
   can produce configuration failures that are difficult to diagnose in
   automated pipelines.

The intended users are ecologists, conservation biologists, and biogeographers
who employ Maxent in reproducible, script-based workflows. `maxentcpp`
addresses these barriers by providing an in-memory, Java-free implementation
that installs as a standard R package with compiled C++ code.

The development version can be installed from GitHub using
`remotes::install_github("alrobles/maxentcpp")`.

# Relationship to existing software

`maxentcpp` fits within a broader modernization trend in ecological software.
Circuitscape.jl [@Hall2021] illustrates how established ecological algorithms
can be reimplemented in a modern high-performance language while retaining
scientific comparability. In R, the spatial ecosystem has moved from older
`sp`/`raster`-centered workflows toward `sf` [@Pebesma2018] and `terra`
[@Hijmans2024], and SDM packages such as `ENMeval` [@Kass2021] and `biomod2`
[@Thuiller2009] have followed this transition.

Within Maxent workflows specifically, existing tools address different parts
of the problem:

- **Java Maxent** [@Phillips2017]: the original implementation and reference
  standard.
- **dismo** [@Hijmans2023] and **kuenm** [@Cobos2019]: R wrappers around Java
  Maxent, inheriting the JDK dependency.
- **maxnet** [@Phillips2024maxnet]: a Java-free R implementation using
  `glmnet` coordinate descent---preserving the same statistical model but
  employing a different optimization algorithm that can produce numerically
  different results, particularly for threshold and hinge features.
- **rmaxent** [@Baumgartner2017]: provides Java-free projection of previously
  fitted Maxent models from `.lambdas` files, but does not replace Java for
  model fitting.
- **ENMeval** [@Kass2021] and **SDMtune** [@Vignali2020]: tuning and
  evaluation frameworks that wrap Maxent-family engines.
- **MIAmaxent** [@Vollering2019]: a related maximum-entropy modeling approach
  using infinitely weighted logistic regression with subset selection rather
  than the Java Maxent optimizer.

`maxentcpp`'s distinguishing contribution is its native C++17 port of the Java
Maxent `density.Sequential` training optimizer, together with per-iteration
validation against the Java implementation via the companion package
`maxentcppCompTest` [@maxentcppCompTest2025].

# Implementation and validation

The package is structured in three layers:

1. **C++17 core** (`src/cpp/include/maxent/`): header-only library implementing
   the `density.Sequential` coordinate-ascent optimizer with $\ell_1$-penalized
   Newton steps, feature transformations, and density normalization. Uses Eigen
   [@Guennebaud2010] for linear algebra.
2. **Rcpp bridge** (`src/`): exposes the C++ engine to R with zero-copy
   semantics where possible.
3. **R interface** (`R/`): user-facing functions for model fitting, spatial
   projection via `terra`, evaluation (AUC, permutation importance), response
   curves, MESS maps, and clamping.

Numerical fidelity is validated by the companion package `maxentcppCompTest`
[@maxentcppCompTest2025], which runs both the C++ optimizer and the original
Java `density.Sequential` on shared test fixtures and compares per-iteration
trajectories of loss, entropy, and $\lambda$ vectors. On controlled fixtures,
agreement reaches $< 10^{-14}$ relative error for symmetric inputs and
$< 10^{-6}$ on $\|\Delta\lambda\|_\infty$ for asymmetric inputs at every
iteration checkpoint.

The package includes an automated test suite run via continuous integration on
Ubuntu, macOS, and Windows with GCC, Clang, and MSVC. A companion vignette
demonstrates ecological equivalence with `maxnet` on a virtual species
(Schoener's $D = 0.93$, Spearman $\rho = 0.95$).

Streaming raster projection processes large rasters in configurable blocks to
bound peak memory, producing results identical to full-memory evaluation under
IEEE 754 double precision with deterministic summation order.

# Limitations

`maxentcpp` currently provides a numerically faithful implementation of the
core continuous-predictor Maxent fitting and projection workflow. It is not a
complete replacement for all Java Maxent functionality. The following features
are not implemented in the current release:

- Categorical predictors
- Replicate runs and cross-validation
- Jackknife variable importance
- Missing-data interpolation
- SWD-to-raster projection mode

Users requiring these features should use Java Maxent via `dismo`, `maxnet`,
or related tools as appropriate. The feature scope is documented in the package
README and pkgdown site.

# AI usage disclosure

Generative AI tools were used during the development of `maxentcpp` and
preparation of this manuscript. GitHub Copilot, Devin, and Claude were used for
code scaffolding, refactoring, test scaffolding, documentation drafting,
continuous-integration configuration, CRAN-compliance review, and manuscript
editing. AI-assisted outputs were reviewed, edited, tested, and validated by
the human author. The human author made the core design decisions, including
the algorithm translation strategy, API design, validation strategy, numerical
fidelity requirements, and release architecture. The author is responsible for
the accuracy, originality, licensing, and maintainability of all submitted
code, documentation, tests, and manuscript text.

# Acknowledgements

The author thanks Steven J. Phillips and collaborators for releasing Maxent as
open-source software, and the Rcpp and RcppEigen developers for enabling
high-performance R/C++ integration.

# References
