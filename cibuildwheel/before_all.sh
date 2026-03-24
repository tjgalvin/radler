#!/bin/bash -eux
#
# This script should be called by `cibuildwheel` in the `before-all` stage.
#
# It installs system packages needed for building Casacore, which in turn is
# needed for building the wheels.

# Set working directory, can be overriden by the user
WORKDIR=${WORKDIR:-/tmp}

# Install system packages needed for building the wheels
function install_system_packages
{
  yum install -y \
    boost1.78-devel \
    cfitsio-devel \
    fftw-devel \
    flex \
    gsl-devel \
    hdf5-devel \
    openblas-devel \
    readline-devel \
    wcslib-devel
}

# Download and build CasaCore using a custom namespace "radler::casacore" to
# avoid conflicts with other CasaCore installations. The environment variable
# CASACORE_VERSION must contain the version number as a dot-separated string,
# e.g. "3.7.1".
function download_and_build_casacore
{
  echo -e "\n==> Downloading and unpacking Casacore ${CASACORE_VERSION} ...\n"
  url="https://github.com/casacore/casacore/archive/refs/tags/v${CASACORE_VERSION}.tar.gz"
  curl -fsSLo - "${url}" | tar -C "${WORKDIR}" -xzf -

  CASACORE_DATA=/usr/local/share/casacore/data
  mkdir -p ${CASACORE_DATA}
  url="https://www.astron.nl/iers/WSRT_Measures.ztar"
  curl -fsSLo - "${url}" | tar -C ${CASACORE_DATA} -xzf -

  echo -e "\n==> Building and installing Casacore ${CASACORE_VERSION} ...\n"
  mkdir -p "${WORKDIR}/casacore-build"
  cd "${WORKDIR}/casacore-build"
  cmake \
    -DCMAKE_CXX_FLAGS="-Dcasacore=radler::casacore" \
    -DBUILD_PYTHON=OFF \
    -DBUILD_PYTHON3=OFF \
    -DBUILD_TESTING=OFF \
    "${WORKDIR}/casacore-${CASACORE_VERSION}"
  make --jobs=`nproc` install
}

set -o pipefail
install_system_packages
download_and_build_casacore
