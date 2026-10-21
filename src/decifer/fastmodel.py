"""
fastmodel.py

Vectorized form of the DeCiFer likelihood.

For a fixed genotype tree (Config) and sample, the VAF of an SNV is an affine
function of the cluster DCF/CCF c (Satas et al. 2021, Eq. 28):

    lam = (c - off) / mu
    v   = (1/F) * (mu * lam + msum)

and c is feasible when it lies within the tree's cell-fraction bounds and lam
lies within [0, lam_limit] (both with a THRESHOLD tolerance). All coefficients
are precomputed here once, for every (mutation, genotype tree) pair ("row")
and every sample, so likelihoods can be evaluated with numpy instead of
per-mutation Python calls to Config.cf_to_v.

The floating-point operations are kept in exactly the same order as
config.Config.cf_to_v and new_coordinate_ascent.compute_pdfs, so the values
are bit-identical to the scalar implementation (see test/test_fastmodel.py).

This module deliberately does not import config/process_input, which parse the
command line at import time.

author: Matt Myers
date: 2026-10-20
"""

import numpy as np
from scipy.special import betaln, gammaln, xlog1py, xlogy

# These must equal config.THRESHOLD and the constants historically defined in
# new_coordinate_ascent.py (checked in test/test_fastmodel.py)
THRESHOLD = 0.05
EPSILON = -1e40
SEQERROR = 1e-40

# maximum number of elements in the temporary arrays of objective_grid
MAX_GRID_ELEMENTS = 250000


class _Params:
    """Coefficients of a set of rows for one sample, ready for broadcasting."""

    __slots__ = ('lbT', 'ubT', 'off', 'mu', 'mu_zero', 'msum', 'invF', 'validF', 'lamlimT',
                 'am1', 'bm1', 'lbeta', 'KS', 'NmK', 'CS', 'SS')


class Model:
    """
    Arrays describing every (mutation, genotype tree) pair of a list of Mutations.

    Rows are ordered by mutation and, within a mutation, by the order of
    mutation.configs; the rows of mutation i are row_start[i]:row_start[i+1].
    Coefficient arrays have shape (num_samples, num_rows).
    """

    def __init__(self, mutations, purity, bb=None):
        self.num_mutations = len(mutations)
        self.num_samples = S = len(mutations[0].a)
        assert all(len(m.a) == S and len(m.d) == S for m in mutations), 'All mutations must have read counts for every sample'
        counts = np.array([len(m.configs) for m in mutations], dtype=np.int64)
        assert np.all(counts > 0), 'Every mutation must have at least one genotype tree'
        self.row_start = np.concatenate(([0], np.cumsum(counts))).astype(np.int64)
        self.row_mut = np.repeat(np.arange(self.num_mutations), counts)
        self.num_rows = R = int(self.row_start[-1])
        self.index_of = {m.index: i for i, m in enumerate(mutations)}
        self.labels = [m.label for m in mutations]
        # log-likelihood columns of the last call of config_cluster_nll, keyed by (sample, cell fraction)
        self._columns = {}

        # Coefficients, computed with the same expressions as config.Config
        lb = np.empty((S, R))
        ub = np.empty((S, R))
        mu = np.empty((S, R))
        msum = np.empty((S, R))
        F = np.empty((S, R))
        lamlim = np.empty((S, R))
        r = 0
        for m in mutations:
            for config in m.configs:
                for s in range(S):
                    # cf_bounds = (M, M + mu), where M is also the offset subtracted in
                    # d_to_lam (DCF) / c_to_lam (CCF): same list of states, same summation order
                    lb[s, r], ub[s, r] = config.cf_bounds(s)
                    mu[s, r] = config.cn_props[config.mut_state][s]
                    msum[s, r] = sum([config.cn_props[st[:2]][s] * st[2] for st in config.other_states])
                    F[s, r] = config.F(s)
                    lamlim[s, r] = purity[s] if config.mut_state == tuple((1, 1)) else 1.0
                r += 1
        self.lbT = lb - THRESHOLD
        self.ubT = ub + THRESHOLD
        self.off = lb
        self.mu = mu
        self.mu_zero = (mu == 0)
        self.msum = msum
        # The scalar code divides by F only for feasible c and would raise ZeroDivisionError;
        # here such rows are simply infeasible
        self.validF = (F != 0)
        with np.errstate(divide='ignore'):
            self.invF = np.where(self.validF, 1. / np.where(self.validF, F, 1.0), np.nan)
        self.lamlimT = lamlim + THRESHOLD

        # Read counts per row, as integer arrays like those built in compute_pdfs
        A = np.array([m.a for m in mutations], dtype=np.int64).T[:, self.row_mut]   # variant reads
        D = np.array([m.d for m in mutations], dtype=np.int64).T[:, self.row_mut]   # total reads
        self.bb = None if bb is None else np.array([bb[s] for s in range(S)], dtype=float)
        if self.bb is None:
            # beta.logpdf(v, a + 1, d - a + 1) = xlog1py(b - 1, -v) + xlogy(a - 1, v) - betaln(a, b)
            AS, BS = A + 1, (D - A) + 1
            self.am1 = AS - 1.0
            self.bm1 = BS - 1.0
            self.lbeta = betaln(AS, BS)
        else:
            KS, BS = A, D - A
            NS = KS + BS
            self.KS = KS.astype(float)
            self.NmK = (NS - KS).astype(float)
            self.CS = gammaln(NS + 1) - gammaln(NS - KS + 1) - gammaln(KS + 1)

    def params(self, s, rows=None, column=False):
        """Coefficients of `rows` (default: all rows) in sample s; column=True shapes them (n, 1)."""
        sel = slice(None) if rows is None else rows
        shape = (-1, 1) if column else (-1,)
        p = _Params()
        for name in ('lbT', 'ubT', 'off', 'mu', 'mu_zero', 'msum', 'invF', 'validF', 'lamlimT'):
            setattr(p, name, getattr(self, name)[s, sel].reshape(shape))
        if self.bb is None:
            p.am1, p.bm1, p.lbeta = (getattr(self, n)[s, sel].reshape(shape) for n in ('am1', 'bm1', 'lbeta'))
        else:
            p.KS, p.NmK, p.CS = (getattr(self, n)[s, sel].reshape(shape) for n in ('KS', 'NmK', 'CS'))
            p.SS = self.bb[s]
        return p

    def assigned_rows(self, mutations):
        """Rows of the genotype trees currently assigned to `mutations`."""
        return np.array([self.row_start[self.index_of[m.index]] + m.configs.index(m.assigned_config) for m in mutations], dtype=np.int64)

    def objective_fn(self, s, rows):
        """
        Returns f(c), equal to new_coordinate_ascent.objective(c, muti, s, bb) where muti are the
        mutations (with their assigned genotype trees) corresponding to `rows`.
        """
        p = self.params(s, rows)
        bb = self.bb is not None

        def f(c):
            v, ok = vaf(c, p)
            if not ok.all():
                return -EPSILON
            return -np.sum(loglik(v, p, bb))
        return f

    def objective_grid(self, cs, s, rows):
        """Values of objective_fn(s, rows) at every point of the 1-D array cs."""
        p = self.params(s, rows)
        bb = self.bb is not None
        cs = np.asarray(cs, dtype=float)
        out = np.full(len(cs), -EPSILON)
        # The objective is -EPSILON unless c is feasible for every row. Row r is feasible only for
        # c in [max(lbT, off - THRESHOLD*mu), min(ubT, off + lamlimT*mu)] (the lam bounds apply when
        # mu > 0), so only points of cs within the intersection of these intervals (with a margin far
        # larger than rounding errors) are evaluated, with the exact kernels.
        with np.errstate(invalid='ignore'):
            lo = np.where(p.mu_zero, p.lbT, np.maximum(p.lbT, p.off - THRESHOLD * p.mu))
            hi = np.where(p.mu_zero, p.ubT, np.minimum(p.ubT, p.off + p.lamlimT * p.mu))
        if not p.validF.all():
            return out
        candidates = np.flatnonzero((cs >= lo.max() - 1e-9) & (cs <= hi.min() + 1e-9))
        step = max(1, MAX_GRID_ELEMENTS // max(1, len(rows)))
        for start in range(0, len(candidates), step):
            idx = candidates[start:start + step]
            v, ok = vaf(cs[idx, None], p)
            feasible = ok.all(axis=1)
            if feasible.any():
                # the log-likelihood is only needed where every row is feasible
                with np.errstate(invalid='ignore'):
                    out[idx[feasible]] = -np.sum(loglik(v[feasible], p, bb), axis=1)
        return out

    def config_cluster_nll(self, C):
        """
        (num_rows, K) array of -sum over samples of the log-likelihood of each (mutation, genotype
        tree) row for each cluster with cell fractions C[sample][cluster]; inf when infeasible.
        """
        bb = self.bb is not None
        tot = None
        # The log-likelihood column of a cluster depends only on the sample and the cluster's cell
        # fraction: columns are computed once per distinct value and kept for the next call, since the
        # absent/truncal/sample-specific clusters are fixed and other centers are often unchanged.
        columns = {}
        for s in range(self.num_samples):
            p = None
            cols = []
            for c in C[s]:
                key = (s, float(c))
                col = columns.get(key)
                if col is None:
                    col = self._columns.get(key)
                if col is None:
                    if p is None:
                        p = self.params(s)
                    v, ok = vaf(c, p)
                    with np.errstate(invalid='ignore'):
                        col = np.where(ok, loglik(v, p, bb), -np.inf)
                columns[key] = col
                cols.append(col)
            ll = np.stack(cols, axis=1)
            # accumulate samples sequentially, as np.sum(axis=0) over a (samples, ...) stack does
            if tot is None:
                tot = ll
            else:
                tot += ll
        self._columns = columns
        return -tot

    def mutation_cluster_nll(self, C):
        """(num_mutations, K) array: best (minimum) value of config_cluster_nll over each mutation's trees."""
        return np.minimum.reduceat(self.config_cluster_nll(C), self.row_start[:-1], axis=0)

    def best_assignment(self, C):
        """
        For every mutation, the (genotype tree, cluster) pair minimizing config_cluster_nll, taking the
        first minimum in (tree, cluster) row-major order like numpy.argmin.

        Returns (nll, best_values, best_rows, best_clusters).
        """
        nll = self.config_cluster_nll(C)
        K = nll.shape[1]
        best_values = np.minimum.reduceat(nll.min(axis=1), self.row_start[:-1])
        flat = np.flatnonzero(nll == best_values[self.row_mut][:, None])
        rows = flat // K
        first = np.flatnonzero(np.diff(self.row_mut[rows], prepend=-1))
        assert len(first) == self.num_mutations, 'Undefined (NaN) likelihood found while assigning mutations: {}'.format(C)
        return nll, best_values, rows[first], flat[first] % K


def vaf(c, p):
    """
    Vectorized config.Config.cf_to_v: returns (v, ok), where v is the VAF clipped to [0, 1] and ok
    is False where cf_to_v would return False. c broadcasts against the arrays of p.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = np.where(p.mu_zero, 0.0, (c - p.off) / p.mu)
        v = p.invF * (p.mu * lam + p.msum)
    ok = (p.lbT <= c) & (c <= p.ubT) & (lam > -THRESHOLD) & (lam < p.lamlimT) & p.validF
    return np.clip(v, 0.0, 1.0), ok


def loglik(v, p, bb):
    """Vectorized new_coordinate_ascent.compute_pdfs for feasible VAFs v (binomial, or beta-binomial if bb)."""
    x = np.minimum(np.maximum(v, SEQERROR), 1.0 - SEQERROR)
    if not bb:
        ll = xlog1py(p.bm1, -x) + xlogy(p.am1, x)
        ll -= p.lbeta
    else:
        SS = p.SS
        minus = 1 - x
        ll = p.CS + betaln(p.KS + SS * x, p.NmK + SS * minus) - betaln(SS * x, SS * minus)
    return np.maximum(ll, EPSILON)
