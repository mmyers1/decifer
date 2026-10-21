"""
decifer.py

author: Simone Zaccaria
date: 2020-05-21
"""

import sys, os
import warnings
import datetime
import traceback
import multiprocessing as mp
import random as rand
from collections import defaultdict

from copy import copy
from multiprocessing import Lock, Value, Pool, Manager
import numpy as np
import math
from bisect import bisect_left

# decifer
from decifer.parse_args import args
from decifer.fileio import write_results, write_results_CIs, read_in_state_trees, write_model_selection_results
from decifer.new_coordinate_ascent import coordinate_descent, objective, compute_pdfs
from decifer.fastmodel import Model
from decifer.mutation import create_mutations
from decifer.process_input import PURITY, MUTATION_DF
from decifer.progress_bar import ProgressBar
from decifer.generator import fit_betabinom

from sklearn.metrics import silhouette_score

def main():
    sys.stderr.write('\n'.join(['Arguments:'] + ['\t{} : {}'.format(a, args[a]) for a in args]) + '\n')

    # create dictionary of sample indices and labels for printing later
    sample_ids = { int(i[0]) : i[1] for i in zip(MUTATION_DF['#sample_index'].unique(), MUTATION_DF['sample_label'].unique()) }
    for i in sample_ids:
        print(i, sample_ids[i])
    num_samples =  len(sample_ids)

    if args['mink'] < 2 + num_samples:
        args['mink'] = 2 + num_samples
        sys.stderr.write('## The minimum number of clusters has been increased to {} to account for fixed clusters!\n'.format(args['mink']))
    if args['maxk'] < args['mink']:
        args['maxk'] = args['mink']
        sys.stderr.write('## The maximum number of clusters has been increased to {} to be higher than the minimum!\n'.format(args['maxk']))
    # state_trees dict: keys are potential CN observed CN states (x,y), e.g. ((1, 0), (1, 1), (0, 0))
    # values are lists of lists, where each list contains all possible genotypes (x,y,m), e.g.
    # [(1, 0, 0), (0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 0, 0)]
    # these lists are always even, and each pair (i,i+1) is an edge in the state tree
    state_trees = read_in_state_trees(args['statetrees'])

    # store info from MUTATION_DF pd.DataFrame in mutations list, containing Mutation objects as elements
    # infer purity from SNV data, but overwrite purity dict if purity file provided (next line)
    mutations = create_mutations(MUTATION_DF, state_trees, not args['ccf'])

    if args['betabinomial']:
        betabinom_out = fit_betabinom(args)
        # invert sample_ids to get index for each sample name, according to input file
        get_index = {v : k for k,v in sample_ids.items()}
        betabinomial = {get_index[sam] : betabinom_out[sam][0] for sam in betabinom_out.keys()}

    # precompute the likelihood coefficients of every mutation and genotype tree (see fastmodel.py)
    model = Model(mutations, PURITY, betabinomial if args['betabinomial'] else None)

    """
    # This is to print state trees and likelihoods for truncal cluster for each mutation, for debugging
    sample = 0 
    for mut in mutations:
        print(mut.label)
        for config,tree in zip(mut.configs, mut.trees):
            tree_toPrint = ';'.join(map(lambda p : '({},{},{})->({},{},{})'.format(*(p[0]+p[1])), zip(tree[:-1:2], tree[1::2])))
            DCF_to_VAF_toPrint = config.cf_to_v(PURITY[sample], sample)


            if args['betabinomial'] is None:
                form = (lambda mut : (config.cf_to_v(PURITY[sample], sample), mut.a[sample] + 1, (mut.d[sample] - mut.a[sample]) + 1))
            else:
                form = (lambda mut : (config.cf_to_v(PURITY[sample], sample), mut.a[sample], mut.d[sample] - mut.a[sample], betabinomial[sample]))
            pdf = -np.sum(compute_pdfs(*zip(*[form(mut)])))
            print(tree_toPrint, DCF_to_VAF_toPrint, pdf)
    sys.exit("Finished printing!!")
    """

    if args['record']:
        manager = mp.Manager()
        record = manager.list()
    if not args['iterative']:
        print("Using binary-search model selection")
        run_coordinator_binary(mutations, model, sample_ids, num_samples, PURITY, args,
                               record if args['record'] else None,
                               betabinomial if args['betabinomial'] else None)
    else:
        print("Using iterative model selection")
        run_coordinator_iterative(mutations, model, sample_ids, num_samples, PURITY, args,
                                  record if args['record'] else None,
                                  betabinomial if args['betabinomial'] else None)
    if args['record']:
        with open('record.log.tsv', 'w') as o:
            o.write('#NUM_CLUSTERS\tRESTART\tSEED\tITERATION\tOBJECTIVE\n')
            for r in record:
                o.write('{}\t{}\t{}\t{}\t{}\n'.format(r[0], r[1], r[2], r[3], r[4]))

def run_coordinator_iterative(mutations, model, sample_ids, num_samples, PURITY, args, record, betabinomial):
    mink, maxk, maxit, prefix, restarts, ubleft, J = unpck(args)
    jobs = [(x, k, np.random.randint(low=0, high=2**10)) for x in range(restarts) for k in range(mink, maxk+1)]
    results = run_descents(jobs, model, num_samples, maxit, record, PURITY, J, args['debug'], args['convergence'])
    # best[cluster number k] = min objective across runs/restarts
    best = {k : min(list(filter(lambda j : j[1] == k, jobs)), key=(lambda j : results[j][4])) for k in range(mink, maxk+1)}

    # Select the best number of clusters given the range of clusters from mink to maxk
    objs = {k: results[best[k]][4] for k in best}
    if args['silhouette']:
        selected, silhouette_scores = silhouette_model_selection({k : unpack_result(mutations, results[best[k]])[1] for k in best}, mink, maxk)
        write_model_selection_results(mink, maxk, objs, silhouette_scores, selected, prefix, args['silhouette'])
    else:
        selected, objs, elbow = elbow_criteria_model_selection(objs, best, mink, maxk, ubleft)
        write_model_selection_results( mink, maxk, objs, elbow, selected, prefix, args['silhouette'])

    if args['printallk']:
        k_to_print = [ k for k in range(mink, maxk+1) ]
    else:
        k_to_print = [ selected ]
    for k in k_to_print:
        finalize_k(mutations, model, results[best[k]], sample_ids, num_samples, PURITY, args, betabinomial, k)

    """
    # FOR TESTING
    #print_PDF(set(clus), bmut, num_samples, args['betabinomial'], C)
    #print_feasibleVAFs(set(clus), bmut, num_samples, args['betabinomial'], C)
    with open("pdfs.txt", 'w') as f:
        for c in set(clus):
            for s in range(num_samples):
                i = str(c) + "_" + str(s) + " "
                f.write(i)
                f.write(" ".join( list(map(str, PDFs[s][c]))))
                f.write("\n")
    with open("max_dcfs.txt", 'w') as f:
        for c in set(clus):
            for s in range(num_samples):
                print c, s, C[s][c], CIs[s][c][0], CIs[s][c][1]
                f.write(" ".join( list(map(str, [c, s, C[s][c], CIs[s][c][0], CIs[s][c][1]] ))))
                f.write("\n")
    """

def finalize_k(mutations, model, result, sample_ids, num_samples, PURITY, args, betabinomial, k):
    # filters poorly fit SNVs, computes cluster CIs and writes the results of the run_descent result for k clusters
    prefix, J = args['output'], args['J']
    C, bmut, clus, conf, objs = unpack_result(mutations, result)
    # C is list of lists; rows are samples, columns are cluster IDs, values are CCFs

    bmut_HQ, bmut_LQ = filter_poorly_fit_SNVs(C, bmut, args['vafdevfilter'])

    #CIs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    CIs, PDFs = compute_CIs_mp(set(clus), bmut_HQ, num_samples, model, J, C, args['debug'], args['conservativeCIs'])

    write_results_CIs(prefix, num_samples, clus, sample_ids, CIs, args['printallk'], k)
    write_results(prefix, C, CIs, conf, bmut_HQ, PURITY, betabinomial, 'CCF' if args['ccf'] else 'DCF', args['printallk'], k, model)
    # write low quality mutations to a separate file
    if len(bmut_LQ) > 0:
        write_results(prefix + "_Outliers", C, CIs, conf, bmut_LQ, PURITY, betabinomial, 'CCF' if args['ccf'] else 'DCF', args['printallk'], k, model)


def silhouette_model_selection(bmut, mink, maxk):
    # bmut[k] = assigned mutations of the best restart for k clusters
    silhouette_scores = {}
    for k in range(mink, maxk + 1):
        mutations = bmut[k]
        X, labels = [], []
        for mut in mutations:
            # don't use mutations that couldn't be assigned to any cluster
            if mut.assigned_cluster > 0:
                labels.append( mut.assigned_cluster )
                var = mut.a
                tot = mut.d
                vaf = [float(v) / t if t > 0 else 0.5 for v, t in zip(var, tot)]
                estC = [mut.assigned_config.v_to_cf(vaf[sam], sam, truncate = False)/PURITY[sam] for sam in range(len(mut.a))]
                X.append(estC)
        silhouette_scores[k] = silhouette_score(np.array(X), np.array(labels))
    selected = max(range(mink, maxk+1), key=(lambda k: silhouette_scores[k]))
    return selected, silhouette_scores

def elbow_criteria_model_selection(objs, best, mink, maxk, ubleft):
    for k in range(mink + 1, maxk + 1):
        if objs[k - 1] < objs[k]:
            best[k] = best[k - 1]
            objs[k] = objs[k - 1]
    chk = (lambda v: v if v != 0.0 else 0.01)
    left = (lambda k: min((objs[k - 1] - objs[k]) / abs(chk(objs[k - 1])), ubleft) if k > mink else ubleft)
    right = (lambda k: (objs[k] - objs[k + 1]) / abs(chk(objs[k])))

    elbow = {k: left(k) - right(k) for k in range(mink, maxk)}
    if mink < maxk:
        selected = max(range(mink, maxk), key=(lambda k: elbow[k]))
    else:
        selected = mink
    return selected, objs, elbow

def filter_poorly_fit_SNVs(C, mut, vafdevfilter):
    # Check for poorly fit SNVs within each sample
    # SNV flagged as low quality (mutation.PASS == 0) if it poorly fits in AT LEAST one sample
    for sample in range(len(C)):
        vaf_deviations, stan_devs = defaultdict(list), defaultdict(float)
        # compute VAF deviations for each mutation
        for m in mut:
            vaf_dev = compute_vaf_dev(C, m, sample)
            vaf_deviations[m.assigned_cluster].append( vaf_dev  )
        # compute standard deviations for each cluster
        for cluster in vaf_deviations:
            stan_devs[cluster] = np.std(np.asarray( vaf_deviations[cluster] ))
        # go back and filter mutations based on normalized VAF deviations
        for m in mut:
            # for the sample-specific clusters, standard deviation will be 0 for some samples
            if stan_devs[m.assigned_cluster] > 0:
                vaf_dev_normalized = compute_vaf_dev(C, m, sample)/stan_devs[m.assigned_cluster]
                if abs(vaf_dev_normalized) > vafdevfilter:
                    m.PASS = 0
    # create list of high and low quality SNVs
    muts_HQ = [m for m in mut if m.PASS == 1]
    muts_LQ = [m for m in mut if m.PASS == 0]
    return muts_HQ, muts_LQ

def compute_vaf_dev(C, mut, sample):
    cluster_vaf = mut.assigned_config.cf_to_v(C[sample][mut.assigned_cluster], sample)
    snv_vaf = float(mut.a[sample]) / float(mut.d[sample])
    return cluster_vaf - snv_vaf

def print_feasibleVAFs(cluster_ids, muts, num_samples, bb, C):
    with open("feasibleVAFs.txt", 'w') as f:
        for i in cluster_ids:
            mut = list(filter(lambda m : m.assigned_cluster == i, muts))
            for s in range(0,num_samples):
                lowers = [m.assigned_config.cf_bounds(s)[0] for m in mut]
                uppers = [m.assigned_config.cf_bounds(s)[1] for m in mut]
                f.write(" ".join(list(map(str, [i, s, max(lowers), min(uppers)]))))
                f.write("\n")

def print_PDF(cluster_ids, muts, num_samples, bb, C):
    with open("pdfs.txt", 'w') as f:
        for i in cluster_ids:
            mut = list(filter(lambda m : m.assigned_cluster == i, muts))
            for s in range(0,num_samples):
                max_dcf = C[s][i] # dcf value that maximizes posterior for this sample and cluster ID
                delta = (-1*objective(max_dcf, mut, s, bb))-2
                #prob = (lambda x: math.exp(-1*(x+delta))) # convert neg log to probability
                prob = (lambda x: math.exp(-1*(x))) # convert neg log to probability
                l = []
                for j in np.linspace(0, 1, 1000):
                   l.append(prob(objective(j, mut, s, bb)))
                total = np.sum(l)
                l = [x/total for x in l]
                f.write(" ".join(list(map(str, [i,s] + l))))
                f.write("\n")


def compute_CIs_mp(cluster_ids, muts, num_samples, model, J, C, debug, conservativeCIs):
    CIs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    PDFs = [[()]*len(C[i]) for i in range(len(C))] # list of lists to store CIs, same structure as C
    num_tests = float(len(cluster_ids)*num_samples) # bonferroni correction for multiple hypothesis testing
    #C[s][i] is the putative mode of the pdf

    CI_fun = CI_conservative if conservativeCIs == True else CI
    jobs = [(c, s, num_tests) for c in cluster_ids for s in range(num_samples)]
    if debug:
        # single-thread mode for debugging
        init_CI(muts, model)
        results = map(CI_fun, jobs)
    else:
        with Pool(processes=min(J, len(jobs)), initializer=init_CI, initargs=(muts, model)) as pool:
            results = list(pool.imap_unordered(CI_fun, jobs))
    for clust, samp, lower, upper, pdf in results:
        CIs[samp][clust] = (lower, upper)
        PDFs[samp][clust] = pdf

    return CIs, PDFs


def init_CI(_muts, _model):
    # mutations (with assignments) and the rows of their assigned genotype trees in the model
    global ci_muts, ci_model, ci_rows, ci_clus
    ci_muts = _muts
    ci_model = _model
    ci_rows = _model.assigned_rows(_muts)
    ci_clus = np.array([m.assigned_cluster for m in _muts], dtype=np.int64)

def CI(job):
    """
    Computes CIs for a sample-cluster combination

    There have been two issues in dealing with the "objective" function to characterize the PDF of CCF/DCF values, needed for obtaining CIs.
    1.) converting -log(unnormalized probability) from objective to an unnormalized probability produced prohibitively large numbers (-log numbers 
    are very negative such that e^-x huge), so we rescaled all values based on the most negative -log value observed from many samples from objective
    (previously we tried using cluster centers, but these do not succeed in finding mode especially when PDF is extremely narrow/disjoint from
    infeasible VAF values truncating distribution).
    2.) across the support of the CCF/DCF PDF distribution, some sample-cluster combinations have 0 values everywhere except an extremely narrow range,
    so trying to integrate approximately with functions like scipy.integrate.quad to get a normalization constant produces 0, because all sampled
    points yield 0.

    Thus, we have used a brute force method where we sample num_pts from the objective function, and use these to create and cahracterize the PDF of
    the CCF/DCF distribution.
    """
    c, s, num_tests = job # c is cluster, s is sample
    rows = ci_rows[ci_clus == c]

    num_pts = 10000
    if len(rows) == 0:
        l, u, pdf = "NaN", "Nan", ["Nan"]*num_pts
    else:
        # objective(j, mut, s, bb) for every grid point j, vectorized
        grid = ci_model.objective_grid(np.linspace(0, PURITY[s], num_pts), s, rows).tolist()
        min_log = min(grid)
        delta = (-1*min_log)-2      # constant to make -log(pdf) values less negative
        prob = (lambda x: math.exp(-1*(x+delta)))           # convert -log(pdf) to unnormalized probability
        total = sum([prob(x) for x in grid])                # unnormalized probabilities across support
        pdf = [prob(x)/total for x in grid]                 # normalized probabilities across support
        cdf = np.cumsum(pdf)

        low_ci = 0.025/num_tests                            # divide the desired CI quantile by the number of tests, bonferonni correction
        high_ci = 1 - low_ci

        low_index = take_closest(cdf, low_ci)
        high_index = take_closest(cdf, high_ci)

        l = float(low_index)/num_pts
        u = float(high_index)/num_pts

    return (c, s, l, u, pdf)


def CI_conservative(job):
    """
    Here we use the DCF point values of the mutations assigned to a cluster to compute that cluster's CIs. We have found
    this produces more conservative estimates of CIs. Specifically, we compute the median of the distribution of DCF point
    values and use bootstrap resampling to calculate CIs. To be conservative, we use the minimum and maximum observed
    median across bootstrap replicates, instead of the medians that correspond to e.g. 0.025 and 0.975 quantiles.
    """
    c, s, num_tests = job  # c is cluster, s is sample
    mut = list(filter(lambda m: m.assigned_cluster == c, ci_muts))
    bootstrap_reps = 10000

    # get DCF point values for cluster
    cf_point_vals = []
    for m in mut:
        vaf = m.a[s]/m.d[s]
        cf_point_vals.append( m.assigned_config.v_to_cf(vaf, s, truncate=False)/PURITY[s] )

    median_boot_reps = bootstrap_medians(cf_point_vals, bootstrap_reps)
    # python floats, so that the CI tuples print the same with numpy 2
    low_CI = float(np.quantile(median_boot_reps, 0.0))
    hi_CI = float(np.quantile(median_boot_reps, 1.0))
    if low_CI == 1.0:
        low_CI = 0.999
    if hi_CI == 0.0:
        hi_CI = 0.001

    return (c, s, low_CI, hi_CI, [None])

def bootstrap_medians(values, reps, max_elements=1000000):
    """
    Medians of `reps` bootstrap resamples of `values`, identical to
        [np.median(rand.choices(values, k=len(values))) for i in range(reps)]
    including the draws from (and the final state of) the global `random` generator, but vectorized:
    numpy's legacy MT19937 RandomState produces the same doubles as random.random() from the same state,
    and rand.choices picks values[floor(random() * n)].
    """
    n = len(values)
    if n == 0:
        return [np.median(rand.choices(values, k=0)) for i in range(reps)]
    values = np.asarray(values, dtype=float)
    version, internal, gauss_next = rand.getstate()
    generator = np.random.RandomState()
    generator.set_state(('MT19937', np.array(internal[:-1], dtype=np.uint32), internal[-1]))
    medians = np.empty(reps)
    step = max(1, max_elements // n)
    for start in range(0, reps, step):
        m = min(step, reps - start)
        # row i holds the n draws of replicate start + i, in the order rand.choices makes them
        idx = np.floor(generator.random_sample((m, n)) * float(n)).astype(np.int64)
        medians[start:start + m] = np.median(values[idx], axis=1)
    _, key, pos, _, _ = generator.get_state()
    rand.setstate((version, tuple(int(k) for k in key) + (int(pos),), gauss_next))
    return medians


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos <= 0:
        return 0
    if pos >= len(myList):
        return pos-1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return pos 
    else:
        return pos-1

def run_coordinator_binary(mutations, model, sample_ids, num_samples, PURITY, args, record, betabinomial):
    L, R, maxit, prefix, restarts, ubleft, J = unpck(args)
    # results[V] = best run_descent result for V clusters; its objective is results[V][4]
    results = {}
    obj = (lambda V : results[V][4])
    def evaluate(V):
        if V not in results:
            sys.stderr.write('[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + 'Computing for {} clusters...\n'.format(V))
            results[V] = run(model, num_samples, V, maxit, PURITY, restarts, J, record, args['debug'], args['convergence'])
        sys.stderr.write('[{:%Y-%b-%d %H:%M:%S}]'.format(datetime.datetime.now()) + 'Objective for {} clusters: {}\n'.format(V, obj(V)))
        return

    MAXR = R
    M = L
    evaluate(MAXR)

    while(R - L > 1):
        M = int(math.floor(float(R + L) / 2.0))
        assert M not in {L, R}
        evaluate(M)
        if float(obj(M) - obj(MAXR)) / abs(obj(MAXR)) > ubleft:
            L = M
        else:
            R = M

    evaluate(L)
    evaluate(R)
    selected = L if float(obj(L) - obj(MAXR)) / abs(obj(MAXR)) <= ubleft else R
    sys.stderr.write('Selected number of clusters: {}\n'.format(selected))
    finalize_k(mutations, model, results[selected], sample_ids, num_samples, PURITY, args, betabinomial, selected)


def run(model, num_samples, K, maxit, PURITY, restarts, J, record, debug, convergence):
    # returns the run_descent result with the best objective across restarts for K clusters
    jobs = [(x, K, np.random.randint(low=0, high=2**10)) for x in range(restarts)]
    # in --debug mode the jobs re-seed the random generators of this process: restore them so that the
    # seeds drawn for the next number of clusters are the same as when the jobs run in a pool
    states = (np.random.get_state(), rand.getstate())
    results = run_descents(jobs, model, num_samples, maxit, record, PURITY, J, debug, convergence)
    np.random.set_state(states[0])
    rand.setstate(states[1])
    return results[min(jobs, key=(lambda j : results[j][4]))]


def run_descents(jobs, model, num_samples, maxit, record, purity, J, debug, convergence):
    # runs the coordinate descent of every job (restart, k, seed); returns {job : run_descent(job)}
    initargs = (model, num_samples, maxit, record, purity, convergence)
    if debug:
        # run in single-thread mode for development/debugging
        init_descent(*initargs)
        return {job : run_descent(job) for job in jobs}
    results = {}
    with Pool(processes=min(J, len(jobs)), initializer=init_descent, initargs=initargs) as pool:
        bar = ProgressBar(total=len(jobs), length=30, verbose=False, lock=Lock(), counter=Value('i', 0))
        bar.progress(advance=False, msg="Started")
        for r in pool.imap_unordered(run_descent, jobs):
            bar.progress(advance=True, msg="Completed {} for k={} [Iterations: {}]".format(r[0], r[1], r[3]))
            results[r[:3]] = r
    return results


def unpack_result(mutations, result):
    # returns C, mutations with the assignments of a run_descent result, cluster and genotype tree assignments, objective
    C, clus, cfg = result[5], result[6].tolist(), result[7].tolist()
    # shallow copies are enough: only the assignments and PASS flag of a copy are ever modified
    bmut = []
    for m, cl, cf in zip(mutations, clus, cfg):
        m = copy(m)
        m.assigned_cluster = cl
        m.assigned_config = m.configs[cf]
        bmut.append(m)
    return C, bmut, clus, [m.assigned_config for m in bmut], result[4]


def init_descent(_model, _num_samples, _maxit, _record, _purity, _convergence):
    global gmodel, num_samples, maxit, record, purity, convergence
    gmodel = _model
    num_samples = _num_samples
    maxit = _maxit
    record = _record
    purity = _purity
    convergence = _convergence


def run_descent(job):
    try:
        x, k, seed = job
        rand.seed(seed)
        np.random.seed(seed)
        #with warnings.catch_warnings() as w:
        #   warnings.simplefilter("ignore")
        np.seterr(all='warn')
        C, clus, cfg, obj, it = coordinate_descent(x, seed, gmodel, num_samples, k, maxit, record, PURITY, convergence)
        return job + (it, obj, C, clus, cfg)
    except:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def unpck(args):
    return args['mink'], args['maxk'], args['maxit'], args['output'], args['restarts'], args['elbow'], args['J']





if __name__ == "__main__":
    main()
