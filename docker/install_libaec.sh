# SPDX-License-Identifier: LGPL-3.0-only

# Script to install libaec (requirement for HDF5)
# libaec implements szip compression, so the optional szip filter can be built
# in HDF5.
set -euo pipefail

pushd /tmp

echo "Downloading libaec"
# The URL includes a hash, so it needs to change if the version does
curl -fsSLO https://gitlab.dkrz.de/k202009/libaec/uploads/ea0b7d197a950b0c110da8dfdecbb71f/libaec-${AEC_VERSION}.tar.gz
tar zxf libaec-${AEC_VERSION}.tar.gz

echo "Building & installing libaec"
pushd libaec-${AEC_VERSION}
./configure
make
make install

# Clean up the files from the build
popd
rm -r libaec-${AEC_VERSION} libaec-${AEC_VERSION}.tar.gz
