"""
new_coordinate_ascent.py

author: Simone Zaccaria
date: 2020-05-25
"""

import sys
import warnings

import numpy as np
from scipy.stats import beta
from scipy.special import betaln
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
import time

from decifer.fastmodel import EPSILON, SEQERROR


TOLERANCE = 1e-03


def coordinate_descent(restart, seed, model, num_samples, num_clusters, MAX_IT=5, record=None, purity=None, convergence='assignments'):
    # model is a fastmodel.Model built from the mutations; assignments are kept as integer arrays:
    # clus[i] is the cluster of mutation i and cfg[i] the index of its genotype tree in mutation.configs
    # convergence: 'assignments' stops when the cluster and genotype-tree assignments of all mutations are
    # unchanged between two consecutive iterations (Satas et al., 2021); 'objective' stops when the objective
    # changed less than TOLERANCE over the last iteration (DeCiFer <= v2.1.4)
    assert convergence in {'assignments', 'objective'}, 'Unknown convergence criterion: {}'.format(convergence)

    # initialize cluster centers C; a list of lists, where each list is the cluster centers for a sample:
    # [absent, truncal, sample specific clusters (truncal if sample index, 0 otherwise),
    # variably present clusters with random initializations]
    if num_samples == 1:
        # if analyzing single sample, don't specify sample-specific truncal clusters, since this is already the case
        # with cluster 1 in single-sample mode
        form = (lambda L, sam: [0.0, purity[sam]] + list(L))
        rest = num_clusters - 2
    else:
        form = ( lambda L, sam : [0.0, purity[sam]] + [purity[sam] if s == sam else 0.0 for s in range(num_samples)] + list(L) )
        # rest = number of variably present clusters subject to model selection, fixed clusters include truncal, absent, and sample specific SNVs
        rest = num_clusters - (2 + num_samples)

    assert rest >= 0, 'Number of clusters is too low!'
    C = [ form(np.random.uniform(low=0.0, high=purity[sam], size=rest) if rest > 0 else [], sam) for sam in range(num_samples) ]
    for sam in range(len(C)):
        assert np.all( np.array(C[sam]) <= purity[sam] ), 'out of bounds!'
    # every mutation starts in cluster 0 with its first genotype tree
    clus = np.zeros(model.num_mutations, dtype=np.int64)
    cfg = np.zeros(model.num_mutations, dtype=np.int64)
    V, C_old, V_old = None, None, None
    IT = 0
    objs = []
    objectives = ClusterObjectives(model, purity)

    # objective criterion: the difference between the last and 3rd-to-last objective less than TOLERANCE
    isconverged = (lambda O : False if len(O) < 3 else abs(O[-1] - O[-3]) < TOLERANCE)
    converged = False
    while not converged and IT < MAX_IT:
        clus_old, cfg_old = clus, cfg
        clus, cfg, obj = optimize_assignments(model, C, clus, cfg)
        assert len(objs) == 0 or obj - objs[-1] <= TOLERANCE, 'Objective is not decreasing after assignments: {}'.format(objs + [obj])
        objs.append(obj)
        if convergence == 'assignments' and IT > 0 and np.array_equal(clus, clus_old) and np.array_equal(cfg, cfg_old):
            # unchanged assignments give unchanged centers for all non-empty clusters: C and the
            # assignments are final
            converged = True
            break
        C_old, V_old = C, update_objs(model, C, clus, cfg, num_samples, num_clusters, objectives)
        C, V, obj = optimize_cluster_centers(model, clus, cfg, num_samples, C_old, V_old, num_clusters, purity, objectives)
        assert len(objs) == 0 or obj - objs[-1] <= TOLERANCE, 'Objective is not decreasing after centers: {}'.format(objs + [obj])
        objs.append(obj)
        IT += 1
        if convergence == 'objective':
            converged = isconverged(objs)
    if not (converged and convergence == 'assignments'):
        clus, cfg, obj = optimize_assignments(model, C, clus, cfg)
        assert len(objs) == 0 or obj - objs[-1] <= TOLERANCE, 'Objective is not decreasing after assignments: {}'.format(objs + [obj])
        objs.append(obj)

    if record is not None:
        for x, obj in enumerate(objs):
            record.append((num_clusters, restart, seed, x, obj))
    return C, clus, cfg, objs[-1], IT


def optimize_assignments(model, C, clus, cfg):
    # For each mutation, select the (genotype tree, cluster) pair with the highest likelihood given the
    # cluster centers C; ties are broken by the first pair in (genotype tree, cluster) order
    nll, best, rows, new_clus = model.best_assignment(C)
    if not np.all(best < np.inf):
        i = int(np.argmax(~(best < np.inf)))
        raise AssertionError('A mutation cannot be assigned to any cluster: {},{}\n\t{}\n'.format(list(model.index_of)[i], model.labels[i], [C[sam][1] for sam in range(model.num_samples)]))
    # the objective of the current assignment can only improve
    before = nll[model.row_start[:-1] + cfg, clus]
    assert np.all(best <= before), 'Non scende: {}'.format(C)
    # python sum over mutations, in mutation order
    return new_clus, rows - model.row_start[:-1], sum(best.tolist())


def optimize_cluster_centers(model, clus, cfg, num_samples, C_old, V_old, num_clusters, purity, objectives):
    # For each sample and cluster, find the DCF/CCF cluster center with minimum -log(prob) given the
    # mutations assigned to the cluster, keeping the previous center if it has a lower objective.
    # Cluster 0 (absent) and 1 (truncal) are fixed; clusters 2..num_samples+1 are sample specific and
    # optimized only in their own sample (0 elsewhere); the remaining clusters are optimized in every
    # sample, or re-initialized at random when empty.
    rows = model.row_start[:-1] + cfg
    members = [rows[clus == i] for i in range(num_clusters)]
    C_new, V_new = [], []
    for sam in range(num_samples):
        C_sam, V_sam = [], []
        for i in range(num_clusters):
            nonempty = len(members[i]) > 0
            if i in {0, 1}:
                x = 0.0 if i == 0 else purity[sam]
            elif i < 2 + num_samples:
                x = objectives.minimum(sam, i, members[i]) if (i - 2) == sam and nonempty else 0.0
            else:
                x = objectives.minimum(sam, i, members[i]) if nonempty else np.random.uniform(low=0.0, high=purity[sam])
            v = objectives.value(sam, i, members[i], x) if nonempty else 0.0
            if V_old[sam][i] < v:
                x, v = C_old[sam][i], V_old[sam][i]
            C_sam.append(x)
            V_sam.append(v)
        C_new.append(C_sam)
        V_new.append(V_sam)
    return C_new, V_new, sum(el for r in V_new for el in r)


def objective(ci, muti, sample, bb):
    # scalar reference implementation of fastmodel.Model.objective_fn, on Mutation objects
    if bb is None:
        form = (lambda mut : (mut.assigned_config.cf_to_v(ci, sample), mut.a[sample] + 1, (mut.d[sample] - mut.a[sample]) + 1))
    else:
        form = (lambda mut : (mut.assigned_config.cf_to_v(ci, sample), mut.a[sample], mut.d[sample] - mut.a[sample], bb[sample]))
    pdfs = compute_pdfs(*zip(*[form(m) for m in muti]), check=True)
    return -EPSILON if np.any(np.isneginf(pdfs)) else -np.sum(pdfs)


def update_objs(model, C, clus, cfg, num_samples, num_clusters, objectives):
    rows = model.row_start[:-1] + cfg
    members = [rows[clus == i] for i in range(num_clusters)]
    iobj = (lambda sam, i : objectives.value(sam, i, members[i], C[sam][i]) if len(members[i]) > 0 else 0.0)
    return [[iobj(sam, i) for i in range(num_clusters)] for sam in range(num_samples)]


class ClusterObjectives:
    """
    Objective (fastmodel.Model.objective_fn) of the mutations assigned to each (sample, cluster). The
    objective of a cluster depends only on its members and center, so its values and its minimizer
    (minimize_scalar is deterministic) are memoized and reused while the members are unchanged.
    """

    def __init__(self, model, purity):
        self.model = model
        self.purity = purity
        # (sample, cluster) -> (rows of the members, objective function, {center : value}, [minimizer])
        self.entries = {}

    def _entry(self, sam, i, rows):
        entry = self.entries.get((sam, i))
        if entry is None or not np.array_equal(entry[0], rows):
            entry = self.entries[(sam, i)] = (rows, self.model.objective_fn(sam, rows), {}, [])
        return entry

    def value(self, sam, i, rows, c):
        _, f, values, _ = self._entry(sam, i, rows)
        key = float(c)
        if key not in values:
            values[key] = f(c)
        return values[key]

    def minimum(self, sam, i, rows):
        # cluster center with minimum objective; min() enforces the purity boundary
        _, f, _, argmin = self._entry(sam, i, rows)
        if not argmin:
            argmin.append(min(minimize_scalar(f, method='bounded', bounds=[0, self.purity[sam]], options={'xatol' : TOLERANCE}).x, self.purity[sam]))
        return argmin[0]


def compute_pdfs(_VS, _AS, _BS, _BB=None, check=False):
     M = np.array(list(map(lambda v : False if v is False else True, _VS)))
     if _BB is None:
         VS, AS, BS = map(np.array, [_VS, _AS, _BS])
     else:
         VS, AS, BS, BB = map(np.array, [_VS, _AS, _BS, _BB])
     VS[VS < SEQERROR] = SEQERROR
     VS[VS > (1.0 - SEQERROR)] = 1.0 - SEQERROR
     RS = np.full(M.shape[0], -np.inf)
     if _BB is None:
        """
        x = beta.logpdf(VS[M], AS[M], BS[M])
        y = VS[M]
        print np.ndim(x)
        print str(x)
        print str(y)
        if np.ndim(x) > 1:
            tmp = str(y) + "\n" + str(x)
            sys.exit(tmp)
        """
        RS[M] = beta.logpdf(VS[M], AS[M], BS[M])
     else:
         NS = AS[M] + BS[M]
         KS = AS[M]
         CS = gammaln(NS + 1) - gammaln(NS - KS + 1) - gammaln(KS + 1) #comb(NS, KS)
         NSminusKS = NS - KS
         BAFS = VS[M]
         minusBAFS = 1 - BAFS
         SS = BB[M]
         RS[M] = CS + betaln(KS + SS * BAFS, NSminusKS + SS * minusBAFS) - betaln(SS * BAFS, SS * minusBAFS)
     RS[M & (RS < EPSILON)] = EPSILON
     if check:
         assert -np.inf not in RS[M]
     return RS

