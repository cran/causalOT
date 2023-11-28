# causalOT 1.0

## Breaking changes
* Simplified the arguments to the main `calc_weight()` function
* Removed direct Python dependency
* Changed estimation function `estimate_effect()` arguments

## New features
* Added a new object-oriented set of functions for advanced users. Described in the new vignette, "Object Oriented COT Objects."
* Added `torch` based data structures allowing GPU support
* Added support for `rkeops` methods to allow for GPU use without the large memory footprint (CRAN version only for Unix machines)
* Added quadratic programming based on the `osqp` package where relevant
* Added internal optimal transport optimization removing reliance on `approxOT` package
* Added `vcov`, `confint`, and `coef` methods for the output of `estimate_effect()`
* Added a barycentric projection outcome estimation method via the function `barycentric_projection()` and a corresponding predict method.
* Added a `plot` method for the `causalWeights` and `causalEffects` objects

## Minor improvements and fixes
* Added a `NEWS.md` file to track changes to the package.