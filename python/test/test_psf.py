# SPDX-License-Identifier: LGPL-3.0-only

import radler as rd
import numpy as np
import pytest

# Tests for the VectorUniquePtrImageAccessor classe. This class is created in
# C++ and can't be created directly by Python. They are available in a
# WorkTableEntry so the tests use that class as a provider for the tested
# classes.


def test_vector_unique_ptr_image_accessor_append():
    data = np.zeros((1, 1), dtype=np.float32)
    entry = rd.WorkTableEntry()
    vector = entry.psfs

    assert len(vector) == 0
    vector.append(data)
    assert len(vector) == 1
