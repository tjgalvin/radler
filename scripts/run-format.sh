#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only

#The directory that contains the source files.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external build CMake)

#The patterns of the C++ source files, which clang-format should format.
CXX_SOURCES=(*.cc *.h)

#The patterns of the CMake source files, which cmake-format should format.
CMAKE_SOURCES=(CMakeLists.txt *.cmake)

#End script configuration.

#The common formatting script has further documentation.
source $(dirname "$0")/../external/aocommon/scripts/format.sh
