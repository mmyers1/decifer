"""
bootstrap_medians (decifer/__main__.py) must reproduce the original bootstrap loop of CI_conservative
exactly, including the state of the global random generator afterwards.
"""

import random as rand

import numpy as np
import pytest


@pytest.mark.parametrize("n", [0, 1, 2, 3, 10, 57, 400])
def test_bootstrap_medians_matches_choices_loop(dec, n):
    import decifer.__main__ as dmain
    values = list(np.random.default_rng(n).normal(size=n) ** 3)
    for seed in [1, 17, 2024]:
        rand.seed(seed)
        rand.random()    # arbitrary position in the generator
        state = rand.getstate()
        ref = [np.median(rand.choices(values, k=len(values))) for i in range(1000)]
        ref_state = rand.getstate()
        rand.setstate(state)
        fast = dmain.bootstrap_medians(values, 1000, max_elements=777)   # several chunks
        assert rand.getstate() == ref_state
        assert np.array_equal(np.asarray(ref), np.asarray(fast), equal_nan=True)
        assert float(np.quantile(ref, 0.0)) == float(np.quantile(fast, 0.0)) or n == 0
