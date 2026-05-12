# Release Checklist for maxentcpp

Use this checklist before each release and JOSS submission.

## Pre-release

- [ ] All tests pass: `devtools::test()`
- [ ] `R CMD check --as-cran` passes with 0 errors, 0 warnings
- [ ] pkgdown site builds: `pkgdown::build_site()`
- [ ] NEWS.md updated with version changes
- [ ] DESCRIPTION version bumped
- [ ] CITATION.cff version and date updated
- [ ] All vignettes render without errors
- [ ] README examples run successfully

## JOSS-specific

- [ ] `paper/paper.md` is within 750–1750 words
- [ ] All references in `paper/paper.bib` are cited in the paper
- [ ] JOSS draft PDF builds via GitHub Action (check Actions artifacts)
- [ ] AI usage disclosure is present and accurate
- [ ] No fabricated claims (CRAN status, benchmarks, adoption)
- [ ] Installation instructions match current state (GitHub vs CRAN)
- [ ] Feature scope claims match implemented features
- [ ] Companion validation package (`maxentcppCompTest`) is accessible

## Release

- [ ] Tag release on GitHub (e.g., `v1.0.0`)
- [ ] Archive on Zenodo (triggers automatically if configured)
- [ ] Record Zenodo DOI
- [ ] Update JOSS review issue with version tag and archive DOI

## Post-acceptance

- [ ] Update CITATION.cff with JOSS DOI
- [ ] Update README with JOSS badge and DOI
- [ ] Update pkgdown site

## JOSS review interaction reminder

> JOSS permits AI use with disclosure, but AI must NOT be used for
> conversational interactions between authors and editors/reviewers
> (except for translation). Answer all JOSS review comments personally.
