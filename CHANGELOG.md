# Radler change log

## Version 1.0.0 (February 11, 2026)

### Improvements

- Support Boost 1.89

### Bug fixes

- Fix calling `radler.perform()` in Python from multiple threads.
- Fix dependency on casacore. Radler requires casacore >= 3.6.0.

## Version 0.1.0 (September 19, 2025)

### Improvements

- Improve documentation
- Require C++20 compiler.
- Modernize building Python binary wheel (using scikit-build-core and cibuildwheel).
- Python: Remove `reached_major_threshold` argument from `radler.perform()`
- Optimize multi-scale gradient descent function.
- Finish major iteration after reaching automask threshold.

### Bug fixes

- Remove dependency on C compiler.
- Fix compilation with AVX

## Version 0.0.1 (August 16, 2022)

Initial release
