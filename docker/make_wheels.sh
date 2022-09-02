#!/usr/bin/bash

# SPDX-License-Identifier: LGPL-3.0-only

# Script to make python wheels for several versions
# The working directory matters, script should be run with ./make_wheels.sh
# The parent directory should contain the source code with the version to build.

set -euo pipefail
for py_version in 310 39 38 37 36; do
    pushd ..
    sed -i "3s/wheel.*/wheel${py_version}/" docker/py310_wheel.docker
    docker build -t radler-py${py_version} -f docker/py310_wheel.docker .
    dockerid=$(docker create radler-py${py_version})
    docker cp $dockerid:/output/ output-${py_version}
    docker rm ${dockerid}
    popd
done
