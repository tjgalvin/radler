# Radler change log

## Version 1.2.0 (October 12, 2026)

### Improvements

- A multiscale RM algorithm was added that performs multiscale CLEAN in Faraday-depth space.

### Bug fixes

- Fix Python binary wheel, which was incompatible with AOFlagger.

## Version 1.1.0 (September 4, 2026)

### Improvements

- There are new options to perform polarimetric deconvolution:
  - A new spectral fitting mode was added to fit spectra to RM values during deconvolution. Enabling this by itself will not perform peak finding in RM space.
  - A new algorithm called "RmDeconvolution" was added, which does perform peak finding in RM space. This can not yet be combined with multiscale ; this is added in version 1.2.0.
  - It is now possible to square channels in combination with multiscale and do deep deconvolutions that minimize the total intensity. In older version, this often led to clean not progressing.
- Combine the point clean `fast_sub_minor_loop` and multi-scale
  `use_sub_minor_optimization` settings into one option.
- Ubuntu 26 is now supported, and Ubuntu 22 is no longer tested.

### Bug fixes

- Fix calling `radler.perform()` in Python from multiple threads.
- Fix main log output when doing sub-image logging.

## Version 1.0.0 (February 11, 2026)

### Improvements

- Support Boost 1.89

### Bug fixes

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
