# swash 3.0.1

## Bugfixes
- Corrections in NAMESPACE (required for correct package installation)
- exponential growth(): R0 is set to NA (instead of 0) if the estimated growth rate is negative
- metrics(): Correcting R-squared for zero variance, added na.rm parameter (default: TRUE)

## Other
- Implemented testthat testing in tests/ folder, keeping former swash_test.R in examples/ folder
- Extensions of documentation


# swash 3.0.0

## Breaking changes (Non-backwards compatible)
- Complete change to the nbmatrix() function: creation of an instance of the new nbmatrix class
- Replacement of the old nbstat() function with the nbstat() method of the nbmatrix class
- Rearrangement of the code in modules for better readability and easier maintaining

## New features
- Spatial statistics based on neighborhood matrix: Global Getis-Ord, Global Moran's I,
- Local Getis-Ord Gi*, with all of them being methods of class nbmatrix
- Plotting a map from an nbmatrix object
- Importing geodata (sf) in infpan objects
- Plotting maps of attributes in an infpan object with method plot_map()

## Bugfixes
- Corrections in RD documentations


# swash 2.0.2

## General
- Deprecation warnings with respect to version >=3.0.0

## Bugfixes
- Updating URLs in help texts that are no longer valid


# swash 2.0.1

## Bugfixes
- Extensions and corrections of documentation files


# swash 2.0.0

## Breaking changes (Non-backwards compatible)
- Analyses are conducted via `infpan` objects rather than `sbm` objects
- Former method `plot_regions()` for class `sbm` is replaced by generic method `plot()` for `infpan` objects
- Former methods `growth()` and `growth_initial()` for `sbm` objects are are now methods for the `infpan` class
- Former function `plot_breakpoints()` is replaced by function `breaks_growth()` and the corresponding `plot()` method

## New features 
- Importing infections panel data via `load_infections_paneldata()`, which creates and instance of the new class `infpan`
- Calculation of spread indicators for `infpan` objects such as incidence and effective reproduction number $R_t$
- Function `hawkes_growth()` for parametrization of a Hawkes process equation for infections
- Method `growth_hawkes()` for parametrization of Hawkes process equations based on `infpan` objects
- Additional NLS estimation in function `exponential_growth()`
- Option to add a constant (if values of y equal to zero occur) in functions `logistic_growth()` and `exponential_growth()`, as well as in methods `growth(infpan)` and `growth_initial(infpan)` 

## Bugfixes
- Functions `metrics()`, `binary_metrics()` and `binary_metrics_glm()` now return results invisible
- Check for equal length of input vectors in `logistic_growth()`


# swash 1.3.3

## New features
- Fit metrics for logistic and exponential growth models

## Bugfixes
- Bug in calculation in metrics() fixed
- Check for vectors lengths in logistic_growth() and exponential_growth()
- Deprecation warnings with respect to version >=2.0.0


# swash 1.3.2

## New features
- New option `verbose` in time-consuming calculation functions.

## Bugfixes
- Deprecation warnings with respect to version >=2.0.0.


