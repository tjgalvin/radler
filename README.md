# Radio Astronomical Deconvolution Library

Radler is a library providing functionality for deconvolving astronomical images. Radler evolved as a stand-alone library from the w-stacking clean (WSClean) imager, `https://gitlab.com/aroffringa/wsclean`_.

## Testing
Radler will be more extensively tested in the near-future. Tests for the core functionality - in particular the different `DeconvolutionAlgorithms` can be found in the `cpp/test` directory. Smaller scale unit tests can be found at namespace level (see, e.g., `cpp/math/test`)

Some example scripts of how the C++ interface can be used, are found in the `cpp/demo` directory.

## License
Radler is released under the LGPL version 3.
