#!/bin/bash
# This script runs the commands to update the strings from a container and updates
# the Python docstring files in your repository to match with the C++ code.
# The clang mkdocs software can be a bit tricky to install in user environment, so
# running this in a container avoids requiring the tools locally.
scriptdir=`dirname $0`
radlerdir=`realpath ${scriptdir}/..`
#sudo docker build -t radler-base -f ubuntu_24_04_base ${radlerdir}/..
cmd=$(cat <<-END
git config --global --add safe.directory /radler &&
rm -rf /radler/docstringbuild &&
mkdir -p /radler/docstringbuild &&
cd /radler/docstringbuild &&
cmake ../ -DBUILD_DOCSTRINGS=ON -DBUILD_PYTHON_BINDINGS=True &&
make docstrings &&
cp doc/docstrings/* /radler/python/docstrings/ &&
cd .. &&
rm -rf docstringbuild/
END
)
sudo docker run --mount src=${radlerdir},target=/radler,type=bind radler-base sh -c "${cmd}"
echo Files in radler/python/docstrings have been updated.

