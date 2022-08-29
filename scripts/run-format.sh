#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only

#The directory that contains the source files.
SOURCE_DIR=$(dirname "$0")/..

#Directories that must be excluded from formatting. These paths are
#relative to SOURCE_DIR.
EXCLUDE_DIRS=(external build CMake python/docstrings)

#End script configuration.

#The common formatting script has further documentation.
source $(dirname "$0")/../external/aocommon/scripts/format.sh
