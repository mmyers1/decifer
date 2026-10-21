"""
Exact-equality tests of the vectorized likelihood (decifer/fastmodel.py) against the scalar
implementation it replaces (config.Config.cf_to_v, new_coordinate_ascent.compute_pdfs/objective,
and the original per-mutation assignment step, reproduced below as an oracle).
"""

import copy

import numpy as np
import pytest


BB = {0: 150.0, 1: 200, 2: 87.5}   # beta-binomial precisions; an int on purpose, as from the brute-force fit


def make(dec, dcf_mode=True, bb=None):
    mutations = dec.mutation.create_mutations(dec.df, dec.state_trees, dcf_mode)
    return mutations, dec.fastmodel.Model(mutations, dec.purity, bb)


def random_C(rng, purity, num_samples, num_clusters):
    return [[0.0, purity[s]] + list(rng.uniform(0.0, purity[s], size=num_clusters - 2)) for s in range(num_samples)]


def random_assignment(rng, mutations, num_clusters):
    muts = [copy.copy(m) for m in mutations]
    for m in muts:
        m.assigned_config = m.configs[rng.integers(len(m.configs))]
        m.assigned_cluster = int(rng.integers(num_clusters))
    return muts


def reference_select(dec, m, C, num_samples, num_clusters, bb):
    # original new_coordinate_ascent.optimize_assignments.select (DeCiFer v2.1.4)
    compute_pdfs = dec.nca.compute_pdfs
    combs = [(config, clust) for config in m.configs for clust in range(num_clusters)]
    if bb is None:
        form = (lambda cf, cl, sam : (cf.cf_to_v(C[sam][cl], sam), m.a[sam]+1, (m.d[sam] - m.a[sam]) + 1))
    else:
        form = (lambda cf, cl, sam : (cf.cf_to_v(C[sam][cl], sam), m.a[sam], m.d[sam] - m.a[sam], bb[sam]))
    grow = (lambda sample : compute_pdfs( *zip(*[form(cb[0], cb[1], sample) for cb in combs]) ))
    objs = -np.sum(np.array([grow(sample) for sample in range(num_samples)]), axis=0)
    return objs, objs.argmin()


def test_constants(dec):
    assert dec.fastmodel.THRESHOLD == dec.config.THRESHOLD
    assert dec.nca.EPSILON == dec.fastmodel.EPSILON == -1e40
    assert dec.nca.SEQERROR == dec.fastmodel.SEQERROR == 1e-40


@pytest.mark.parametrize("dcf_mode", [True, False])
def test_vaf_matches_cf_to_v(dec, dcf_mode):
    mutations, model = make(dec, dcf_mode)
    S = model.num_samples
    checked = infeasible = 0
    r = 0
    for m in mutations:
        for config in m.configs:
            for s in range(S):
                lb, ub = config.cf_bounds(s)
                mu = config.cn_props[config.mut_state][s]
                cs = list(np.linspace(0.0, dec.purity[s], 41)) + [lb, ub, lb - 0.05, ub + 0.05, lb - 0.05 * mu, ub + 0.05 * mu]
                cs += [np.nextafter(c, d) for c in cs[41:] for d in (-np.inf, np.inf)]
                v, ok = dec.fastmodel.vaf(np.array(cs), model.params(s, [r]))
                for j, c in enumerate(cs):
                    ref = config.cf_to_v(float(c), s)
                    if ref is False:
                        assert not ok[j], (r, s, c)
                        infeasible += 1
                    else:
                        assert ok[j] and ref == v[j], (r, s, c, ref, v[j])
                    checked += 1
            r += 1
    assert infeasible > 0 and checked > infeasible


@pytest.mark.parametrize("bb", [None, BB])
@pytest.mark.parametrize("dcf_mode", [True, False])
def test_loglik_matches_compute_pdfs(dec, dcf_mode, bb):
    mutations, model = make(dec, dcf_mode, bb)
    rows = [(m, config) for m in mutations for config in m.configs]
    for s in range(model.num_samples):
        for c in np.linspace(0.0, dec.purity[s], 13):
            if bb is None:
                args = zip(*[(config.cf_to_v(c, s), m.a[s] + 1, (m.d[s] - m.a[s]) + 1) for m, config in rows])
            else:
                args = zip(*[(config.cf_to_v(c, s), m.a[s], m.d[s] - m.a[s], bb[s]) for m, config in rows])
            ref = dec.nca.compute_pdfs(*args)
            p = model.params(s)
            v, ok = dec.fastmodel.vaf(c, p)
            fast = np.where(ok, dec.fastmodel.loglik(v, p, bb is not None), -np.inf)
            assert np.array_equal(ref, fast)


@pytest.mark.parametrize("bb", [None, BB])
@pytest.mark.parametrize("dcf_mode", [True, False])
def test_objective_matches_scalar(dec, dcf_mode, bb):
    rng = np.random.default_rng(1)
    mutations, model = make(dec, dcf_mode, bb)
    num_clusters = 6
    for _ in range(3):
        muts = random_assignment(rng, mutations, num_clusters)
        rows = model.assigned_rows(muts)
        clus = np.array([m.assigned_cluster for m in muts])
        for i in range(num_clusters):
            muti = [m for m in muts if m.assigned_cluster == i]
            if not muti:
                continue
            for s in range(model.num_samples):
                cs = np.linspace(0, dec.purity[s], 57)
                ref = [dec.nca.objective(c, muti, s, bb) for c in cs]
                f = model.objective_fn(s, rows[clus == i])
                assert [f(c) for c in cs] == ref
                assert model.objective_grid(cs, s, rows[clus == i]).tolist() == ref


@pytest.mark.parametrize("bb", [None, BB])
@pytest.mark.parametrize("dcf_mode", [True, False])
def test_assignment_matches_reference(dec, dcf_mode, bb):
    rng = np.random.default_rng(2)
    mutations, model = make(dec, dcf_mode, bb)
    S, K = model.num_samples, 7
    for _ in range(4):
        C = random_C(rng, dec.purity, S, K)
        nll, best_values, best_rows, best_clus = model.best_assignment(C)
        per_mut = model.mutation_cluster_nll(C)
        for i, m in enumerate(mutations):
            objs, best = reference_select(dec, m, C, S, K, bb)
            start, end = model.row_start[i], model.row_start[i + 1]
            assert np.array_equal(nll[start:end].ravel(), objs)
            assert best_rows[i] - start == best // K and best_clus[i] == best % K
            assert best_values[i] == objs[best]
            # original fileio.write_results "LHs": minimum over genotype trees for each cluster
            assert np.array_equal(per_mut[i], objs.reshape(-1, K).min(axis=0))


def test_reduction_orders():
    # assumptions behind bit-identical sums: (1) row sums of a C-contiguous 2-D array equal the 1-D
    # (pairwise) sums of each row; (2) np.sum(axis=0) over a stack accumulates the stack sequentially
    rng = np.random.default_rng(3)
    for n in [1, 7, 8, 9, 100, 129, 1000, 8193, 20000]:
        a = rng.normal(size=(5, n)) * 10.0 ** rng.integers(-5, 40, size=(5, n))
        assert np.array_equal(a.sum(axis=1), np.array([np.sum(row) for row in a]))
    for S in [1, 2, 3, 8, 11, 17]:
        stack = rng.normal(size=(S, 50, 7)) * 10.0 ** rng.integers(-5, 40, size=(S, 50, 7))
        acc = stack[0].copy()
        for s in range(1, S):
            acc += stack[s]
        assert np.array_equal(np.sum(stack, axis=0), acc)
