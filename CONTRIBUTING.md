# Contributing to swash

Thank you for your interest in contributing to `swash`.

`swash` is an open-source R package for quantitative analysis in health geography, with a particular focus on the analysis of spatial spreads of infectious diseases based on infections panel data. The package provides methods for analysing epidemic growth, spatial diffusion, and hotspots/clusters, including the Swash-Backwash Model for the Single Epidemic Wave, phenomenological growth models, and spatial statistics.

## Reporting issues

Bug reports, questions, and feature requests are welcome.

Please use the [GitHub Issues](https://github.com/geowieland/swash_official/issues) page to report:

* bugs or unexpected behaviour,
* documentation issues,
* feature requests,
* questions about the use of the package.

When reporting a bug, please provide a minimal reproducible example where possible, including the `swash` version, relevant R and package versions, and information about the data or input that caused the problem.

## Contributing code

Contributions via pull requests are welcome.

Before submitting a pull request:

1. Fork the repository and create a separate branch for your changes.
2. Keep changes focused on a single feature, bug fix, or documentation improvement.
3. Follow the existing code structure and style of the package.
4. Update the documentation and examples when appropriate.
5. Add or update tests where appropriate.
6. Make sure that existing functionality is not unintentionally changed.
7. Check that the package can be built and checked successfully before submitting the pull request.

Pull requests should include a short description of the changes and their motivation.

## Documentation

Improvements to the documentation, examples, and methodological explanations are welcome.

Documentation changes should be consistent with the existing terminology and structure of the project. Examples should, where appropriate, use the existing `infpan` framework and demonstrate reproducible workflows.

## Development

The package can be installed directly from CRAN:

```r
install.packages("swash")
```

For development, clone the repository and install the package locally:

```bash
git clone https://github.com/geowieland/swash_official.git
cd swash_official
```

The package can then be installed from the local repository, for example using `remotes`:

```r
install.packages("remotes")
remotes::install_local()
```

The `tests/` directory contains usage examples and tests for many of the included functions.

When developing new functionality, please consider both the underlying methodological implementation and its integration with the existing `infpan` class and associated `summary()` and `plot()` methods where applicable.

## Scientific contributions

Contributions related to the implementation or extension of epidemiological, statistical, spatial, or health-geographical methods should include appropriate references to the underlying scientific literature where relevant.

Please describe methodological changes clearly so that their scientific purpose, assumptions, and implementation can be reviewed. In particular, contributions involving model-based analyses should clearly document the methodological basis of the implementation.

If you use the `swash` R package in scientific work, please cite the software according to the citation information provided in the repository.

## Code of conduct

Please keep discussions and contributions respectful, constructive, and focused on improving the software.

## License

By contributing to this repository, you agree that your contributions will be licensed under the GPL-3.0 License used by `swash`.