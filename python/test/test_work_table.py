import numpy as np
import pytest

import radler as rd


def test_work_table_entry():
    entry = rd.WorkTableEntry()

    # Check defaults
    assert entry.image_weight == 0.0
    assert entry.central_frequency == 0.0
    assert entry.band_start_frequency == 0.0
    assert entry.band_end_frequency == 0.0
    assert entry.original_channel_index == 0
    assert entry.original_interval_index == 0
    assert entry.mask_channel_index == 0

    # Update defaults
    entry.image_weight = 1.25
    entry.band_start_frequency = 50.0e6
    entry.band_end_frequency = 60.0e6
    entry.original_channel_index = 2
    entry.original_interval_index = 1
    entry.mask_channel_index = 42

    # Check updated values
    assert entry.image_weight == 1.25
    assert entry.central_frequency == (50.0e6 + 60.0e6) / 2.0
    assert entry.band_start_frequency == 50.0e6
    assert entry.band_end_frequency == 60.0e6
    assert entry.original_channel_index == 2
    assert entry.original_interval_index == 1
    assert entry.mask_channel_index == 42


def test_zero_groups():
    """
    Check WorkTable constructor with zero original / deconvolution groups.
    """

    n_original_groups = 0
    n_deconvolution_groups = 0
    work_table = rd.WorkTable([], n_original_groups, n_deconvolution_groups)
    assert work_table.original_groups == [[]]
    assert work_table.deconvolution_groups == [[0]]


def test_negative_original_groups():
    """
    Check WorkTable constructor for negative number
    of original groups.
    """

    n_original_groups = -2
    n_deconvolution_groups = 1
    with pytest.raises(TypeError):
        work_table = rd.WorkTable(
            [], n_original_groups, n_deconvolution_groups
        )


def test_negative_deconvolution_groups():
    """
    Check WorkTable constructor for negative number
    of deconvolution groups.
    """

    n_original_groups = 10
    n_deconvolution_groups = -1
    with pytest.raises(TypeError):
        work_table = rd.WorkTable(
            [], n_original_groups, n_deconvolution_groups
        )


@pytest.mark.parametrize(
    "n_original_groups,n_deconvolution_groups",
    [(4, 12), (12, 4)],
)
def test_multiple_deconvolution_groups(
    n_original_groups, n_deconvolution_groups
):
    """
    Test the WorkTable constructor for (combinations of) multiple deconvolution groups.
    """

    work_table = rd.WorkTable([], n_original_groups, n_deconvolution_groups)

    assert len(work_table.original_groups) == n_original_groups
    assert not any(work_table.original_groups)
    assert len(work_table.deconvolution_groups) == min(
        n_deconvolution_groups, n_original_groups
    )

    if n_original_groups > 1:
        n_sub_groups = n_original_groups // min(
            n_original_groups, n_deconvolution_groups
        )
        for i, sub_table in enumerate(work_table.deconvolution_groups):
            ref_list = list(range(i * n_sub_groups, (i + 1) * n_sub_groups))
            assert sub_table == ref_list


@pytest.mark.parametrize("channel_index_offset", [-3, None, 2])
def test_channel_index_offset(channel_index_offset):
    """
    Test the value of channel_index_offset in the constructor.
    """

    n_original_groups = 2
    n_deconvolution_groups = 2

    if channel_index_offset is not None and channel_index_offset < 0:
        with pytest.raises(TypeError):
            work_table = rd.WorkTable(
                [],
                n_original_groups,
                n_deconvolution_groups,
                channel_index_offset,
            )
    else:
        work_table = (
            rd.WorkTable([], n_original_groups, n_deconvolution_groups)
            if channel_index_offset is None
            else rd.WorkTable(
                [],
                n_original_groups,
                n_deconvolution_groups,
                channel_index_offset,
            )
        )
        assert (
            work_table.channel_index_offset == 0
            if channel_index_offset is None
            else channel_index_offset
        )


def test_psfs_wrong_shape():
    """
    Check that invalid PSF entries throw an exception.
    """

    with pytest.raises(TypeError):
        rd.WorkTable([1], 1, 1)  # Not a 2D array

    with pytest.raises(TypeError):
        rd.WorkTable([1, 1], 1, 1)  # Not a 2D array

    with pytest.raises(TypeError):
        rd.WorkTable([[1, 1], []], 1, 1)

    with pytest.raises(TypeError):
        rd.WorkTable([[1, 1], [1]], 1, 1)

    with pytest.raises(TypeError):
        rd.WorkTable([[1, 1], [1, 1, 1]], 1, 1)


def test_add_entries_wrong_type():
    """
    Check that TypeErrors are thrown in case provided
    images are not np.float32
    """

    entries = [rd.WorkTableEntry(), rd.WorkTableEntry(), rd.WorkTableEntry()]

    psf = np.ones((4, 4), dtype=np.float64)
    residual = np.ones((4, 4), dtype=int)
    model = np.ones((4, 4), complex)

    with pytest.raises(TypeError):
        entries[0].psfs.append(psf)

    with pytest.raises(TypeError):
        entries[1].residual = residual

    with pytest.raises(TypeError):
        entries[2].model = model


def test_add_entries():
    """
    Check the add_entry member function.
    Largely following the add_entries test in test_work_table.cc.

    NOTE: just checking whether the interface behaves correctly.
    """
    work_table = rd.WorkTable([], 3, 1)

    entries = [rd.WorkTableEntry(), rd.WorkTableEntry(), rd.WorkTableEntry()]

    psf = np.ones((4, 4), dtype=np.float32)
    residual = np.ones((4, 4), dtype=np.float32)
    model = np.ones((4, 4), np.float32)

    # Write property
    entries[0].psfs.append(psf)
    entries[1].residual = residual
    entries[2].model = model

    # Read (image) properties back should fail
    with pytest.raises(AttributeError):
        residual_read = entries[1].residual

    with pytest.raises(AttributeError):
        model_read = entries[1].model

    assert len(work_table) == 0
    for i, entry in enumerate(entries):
        entry.image_weight = i
        entry.original_channel_index = i % 2
        entry.band_end_frequency = float(i + 1) * 1e6
        work_table.add_entry(entry)
    assert len(work_table) == len(entries)

    assert len(work_table.original_groups) == len(entries)
    assert len(work_table.original_groups[0]) == 2
    assert len(work_table.original_groups[1]) == 1
    assert len(work_table.original_groups[2]) == 0

    # Check iterator
    for i, entry in enumerate(work_table):
        assert entry.image_weight == i
        assert entry.central_frequency == (0.5 * float(i + 1) * 1e6)


def test_str(capsys):
    """
    Check the WorkTable.__str__().

    Both str(WorkTable) and print(WorkTable) use the __str__ member function
    and are tested to validate they work as intended.
    """

    work_table = rd.WorkTable([[1, 2], [3, 4], [5, 6]], 0, 0)
    assert (
        work_table.__str__()
        == """=== IMAGING TABLE ===
Original groups       1
Deconvolution groups  1
Channel index         0
=== PSFs ===
[x: 1, y: 2]
[x: 3, y: 4]
[x: 5, y: 6]
"""
    )
    assert (
        str(work_table)
        == """=== IMAGING TABLE ===
Original groups       1
Deconvolution groups  1
Channel index         0
=== PSFs ===
[x: 1, y: 2]
[x: 3, y: 4]
[x: 5, y: 6]
"""
    )

    print(work_table)
    captured = capsys.readouterr()
    assert (
        captured.out
        == """=== IMAGING TABLE ===
Original groups       1
Deconvolution groups  1
Channel index         0
=== PSFs ===
[x: 1, y: 2]
[x: 3, y: 4]
[x: 5, y: 6]

"""
    )
