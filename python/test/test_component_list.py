import numpy as np
import pytest

import radler as rd


def test_component_list():
    cl = rd.ComponentList()

    # Check defaults
    assert cl.width == 0
    assert cl.height == 0
    assert cl.n_scales == 0
    assert cl.n_frequencies == 0

    with pytest.raises(IndexError):
        cl.component_count(0)


def test_multiscale_component_list():
    WIDTH = 8
    HEIGHT = 9
    N_SCALES = 3
    N_FREQUENCIES = 7

    cl = rd.ComponentList(WIDTH, HEIGHT, N_SCALES, N_FREQUENCIES)

    assert cl.width == WIDTH
    assert cl.height == HEIGHT
    assert cl.n_scales == N_SCALES
    assert cl.n_frequencies == N_FREQUENCIES

    for scale_num in range(N_SCALES):
        assert cl.component_count(scale_num) == 0

    with pytest.raises(IndexError):
        cl.component_count(N_SCALES)
