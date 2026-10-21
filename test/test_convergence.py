"""
Tests of the stopping rules of new_coordinate_ascent.coordinate_descent.
"""

import numpy as np
import pytest


def run(dec, monkeypatch, convergence, seed, num_clusters=7, maxit=200):
    mutations = dec.mutation.create_mutations(dec.df, dec.state_trees, True)
    model = dec.fastmodel.Model(mutations, dec.purity, None)
    # record the assignments returned by every assignment step
    calls = []
    original = dec.nca.optimize_assignments
    def spy(*args):
        out = original(*args)
        calls.append((out[0].copy(), out[1].copy()))
        return out
    monkeypatch.setattr(dec.nca, 'optimize_assignments', spy)
    np.random.seed(seed)
    result = dec.nca.coordinate_descent(0, seed, model, model.num_samples, num_clusters, maxit, None, dec.purity, convergence)
    return result, calls


same = (lambda x, y : np.array_equal(x[0], y[0]) and np.array_equal(x[1], y[1]))


@pytest.mark.parametrize("seed", range(6))
def test_assignments_rule_stops_at_first_repeat(dec, monkeypatch, seed):
    (C, clus, cfg, obj, it), calls = run(dec, monkeypatch, 'assignments', seed)
    repeats = [t for t in range(1, len(calls)) if same(calls[t], calls[t - 1])]
    # stops right after the first assignment step that repeats the previous one, without a final extra step
    assert repeats == [len(calls) - 1]
    assert it == len(calls) - 1
    assert same((clus, cfg), calls[-1])


def test_assignments_rule_respects_maxit(dec, monkeypatch):
    (C, clus, cfg, obj, it), calls = run(dec, monkeypatch, 'assignments', 0, maxit=1)
    assert it == 1 and len(calls) == 2


def test_unknown_rule(dec, monkeypatch):
    with pytest.raises(AssertionError):
        run(dec, monkeypatch, 'bogus', 0)
