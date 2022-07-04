# SPDX-License-Identifier: LGPL-3.0-only

import radler as rd
import numpy as np
import pytest

# Tests for the VectorPsf and Psf classes. These classes are created in C++ and
# can't be created directly by Python. They are available in a WorkTableEntry
# so the tests use that class as a provider for the tested classes.


def test_vector_psf_append():
    data = np.zeros((1, 1), dtype=np.float32)
    entry = rd.WorkTableEntry()
    vector = entry.psfs

    # 1 argument
    vector.append(data)
    psf = vector[0]
    assert psf.x == 0
    assert psf.y == 0
    with pytest.raises(AttributeError):  # Can't read
        dummy = psf.accessor

    # 3 arguments
    vector.append(42, 99, data)
    psf = vector[1]
    assert psf.x == 42
    assert psf.y == 99
    with pytest.raises(AttributeError):  # Can't read
        dummy = psf.accessor


def test_vector_psf_assignment():
    data = np.zeros((1, 1), dtype=np.float32)
    entry = rd.WorkTableEntry()
    vector = entry.psfs

    vector.append(data)
    psf = vector[0]

    psf.x = 42
    assert psf.x == 42
    assert psf.y == 0

    psf.y = 99
    assert psf.x == 42
    assert psf.y == 99

    psf.accessor = data
    assert psf.x == 42
    assert psf.y == 99


def test_vector_psf_print():
    data = np.zeros((1, 1), dtype=np.float32)
    entry = rd.WorkTableEntry()
    vector = entry.psfs

    assert str(vector) == "[]"

    vector.append(data)
    psf = vector[0]
    assert str(psf) == "[x: 0, y: 0, accessor: set]"
    assert str(vector) == "[[x: 0, y: 0, accessor: set]]"

    vector.append(42, 99, data)
    psf = vector[1]
    assert str(psf) == "[x: 42, y: 99, accessor: set]"
    assert str(vector) == "[[x: 0, y: 0, accessor: set], [x: 42, y: 99, accessor: set]]"


def test_vector_psf_size_and_access_validation():
    data = np.zeros((1, 1), dtype=np.float32)
    entry = rd.WorkTableEntry()
    vector = entry.psfs

    assert len(vector) == 0
    with pytest.raises(IndexError):  # Index out of bounds
        vector[0]

    for i in range(1, 5):
        vector.append(data)
        assert len(vector) == i
        vector[i - 1]
        with pytest.raises(IndexError):  # Index out of bounds
            vector[i]
        with pytest.raises(IndexError):  # Index out of bounds
            vector[i + 1]
