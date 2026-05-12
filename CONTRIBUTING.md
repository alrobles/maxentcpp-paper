# Contributing to maxentcpp

Thank you for your interest in contributing to `maxentcpp`.

## Reporting bugs

Please open an issue at <https://github.com/alrobles/maxentcpp/issues> with:

- A minimal reproducible example
- Output of `sessionInfo()`
- Operating system and R version
- Expected vs. actual behavior

## Feature requests

Open an issue describing the use case and expected behavior. If the feature
relates to a Java Maxent capability not yet implemented, please note the
Java Maxent version and settings involved.

## Running tests

```r
# Install development dependencies
remotes::install_deps(dependencies = TRUE)

# Run the test suite
devtools::test()

# Full R CMD check
devtools::check()
```

## Pull requests

1. Fork the repository and create a feature branch from `main`.
2. Add tests for any new functionality.
3. Run `devtools::check()` and ensure 0 errors, 0 warnings.
4. Update documentation if needed (`devtools::document()`).
5. Submit a pull request with a clear description of the change.

## Code style

- Follow existing conventions in the codebase.
- R code: use roxygen2 for documentation, `testthat` for tests.
- C++ code: C++17 standard, header-only where possible, no `std::cout`.

## Maintainer response

The maintainer aims to respond to issues and PRs within one week. Complex
contributions may take longer for review.
