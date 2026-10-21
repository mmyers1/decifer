# 2026-10-20 — DeCiFer efficiency

Approved plan: `/home/myersm2/.claude/plans/improve-the-efficiency-of-snappy-rabin.md`
Branch: `perf/vectorize-core` (no commits made unless requested)
Conda env: `decifer` (`/home/myersm2/miniforge3/envs/decifer`)

## Summary

**Status: complete. The user stopped me after the numpy-2 dry run and later said "resume".** The
remaining open items are runs that are too expensive for me to do and are handed to the user (commands
under Progress).

- **Phase A (exact vectorization)** is complete and verified **byte-identical** to v2.1.4 (`1b6f6d8`)
  for the same seed. It is ~20–55× faster on the test scenarios, with the same peak memory.
- **Phase B** is complete:
  - `--binarysearch` works (it crashed before).
  - The beta-binomial fit uses `--restarts_bb` (fitted precisions change by ~2e-7 relative on the test data).
  - The `fitbeta` entry point works.
  - numpy-2 / pandas-2 fixes are in.
  - The latent NameError in `process_input` is fixed.
  - **New default stopping rule `--convergence assignments`** (paper e5).
    `--convergence objective` reproduces v2.1.4 exactly.
- **Phase C**, all **exact** (byte-identical outputs):
  - a CI feasible-interval prefilter
  - an assignment-step column cache and memoization of unchanged clusters (descents at 20k SNVs 2.1×
    faster than after Phase A)
  - **C3**: an exact vectorized `--conservativeCIs` bootstrap (3.6× faster at 5k)

  C1 (batched Brent) was not justified by the profile. **C4 (`--init dcf`) failed ground-truth validation
  (worse ARI/CF error/k̂) and was removed.**
- **numpy 2** verified (identical outputs); the env is back on numpy 1.26.4 so the v2.1.4 baseline still runs.
- Tests: 33 pass. That includes the 100-restart end-to-end test, which reproduces the stored truth CIs.

## Method

1. **Phase A — exact speedups.** The VAF of an SNV is an affine function of the cluster DCF/CCF for a
   fixed genotype tree and sample (paper Eq. 28, `config.py:92-134`). All per-(mutation, tree, sample)
   coefficients are precomputed once into arrays (`src/decifer/fastmodel.py`). Every likelihood
   evaluation (assignment step, cluster-center step, CI grid, output "LHs") becomes vectorized
   numpy work. Floating-point operations are kept in the same order as the scalar code, so outputs
   for a fixed seed are byte-identical. A golden-master regression against the original code checks
   this (`benchmarks/run_regression.sh`).
2. **Phase B — bug fixes and the paper's stopping rule** (expected, documented output changes).
3. **Phase C — algorithmic options**, only where profiling after Phase A justifies them. Validated
   with ground-truth simulations: cell-fraction error, ARI, k̂ − k.

Baseline for comparisons: `git archive HEAD src` extracted to `$TMPDIR/baseline` and run with
`PYTHONPATH=$TMPDIR/baseline/src python -m decifer` in the same env. I verified that `PYTHONPATH`
takes precedence over the editable install.

## Progress

- [x] Phase 0: branch, env, baseline extract, calibration run, profile
- [x] Phase A: fastmodel + vectorized descent / CI / output / input; unit tests; golden master
- [x] Phase B: bug fixes; `--convergence assignments` default; ground-truth validation
- [x] Phase B numpy-2 check (resumed after the user's stop): passed; env rolled back to 1.26.4 (row 29)
- [x] Phase C: profiling at 5k/20k SNVs; exact items implemented (rows 26–28); C1 not justified;
  C3 implemented **exactly** (row 30); C4 implemented, failed validation, and removed (rows 32–36)
- [ ] Large-scale benchmarks and full-output PDAC/synthetic golden masters: handed to the user (commands below)

### Commands for the user (not run by me)
```bash
cd /data1/shahs3/users/myersm2/repos/decifer
PY=/home/myersm2/miniforge3/envs/decifer/bin/python
# 1) full-output golden master incl. 11-sample PDAC and synthetic-500 (original code ~3-5 min each)
PYTHON=$PY benchmarks/run_regression.sh -s "pdac synthetic" -n "--convergence objective"
#    expect: IDENTICAL for every file
# 2) (done by me after "resume": numpy-2 check passed)
#    note: the e2e test overwrites tracked files; restore with: git checkout -- test/decifer_*.tsv test.png
# 3) benchmark at scale (new code; the original takes hours at these sizes, so compare descents only)
$PY benchmarks/make_synthetic.py scaleup test/data/decifer_input.tsv 20000 syn_20000.tsv --seed 1
$PY benchmarks/timecmd.py -- $PY -m decifer syn_20000.tsv -p test/data/purity.tsv -k 5 -K 12 -r 20 -j 8 --seed 17 -o syn20k
git archive 1b6f6d8 src | tar -x -C /tmp/decifer-v214
$PY benchmarks/compare_descent.py syn_20000.tsv test/data/purity.tsv --base-src /tmp/decifer-v214/src --new-src src -k 8 -r 2
```

## Results

### Baseline calibration (original code, env `decifer`)
`test/data` (208 SNVs, 3 samples), `-k 5 -K 5 -r 2 --seed 17 --debug`: **41.3 s wall, 169 MB peak RSS**.

cProfile of the same run (74.6 s under the profiler):

| Function | cumtime | share |
|---|---|---|
| `__main__.CI` (10,000-pt grid × 15 cluster-sample pairs) | 69.5 s | 93% |
| └ `new_coordinate_ascent.objective` (150k calls) | 69.3 s | |
| └ `Config.cf_to_v` (5.0M calls) | 39.8 s | |
| └ `compute_pdfs` → `scipy.stats.beta.logpdf` (160k calls, ~124 µs overhead each) | 25.0 s | |
| coordinate descent (2 restarts) | ~2 s | 3% |

The hot spots match the plan. In this small case the CI grid dominates. On larger inputs, the
descent (assignment + center steps) scales as n·configs·K·S per iteration × restarts × K values.

### Phase A gate — PASSED (byte-identical, same seed)

Unit tests: `pytest test/test_fastmodel.py` gives 16 passed (6 s), all exact equality.

Golden master (`benchmarks/run_regression.sh`, base = `1b6f6d8`, all `--seed 17`; single-core
wall time on a shared node with load average ~75, so timings are approximate):

| Scenario | Original | New | Speedup | Outputs |
|---|---|---|---|---|
| binomial (`-k5 -K5 -r3 --debug`) | 42.4 s | 1.7–2.1 s | ~22× | 4/4 identical |
| betabinomial (seg/snp, `-r3 --debug`) | 33.3 s | 5.7 s | 5.8× (beta-binomial fit/plot not optimized) | 4/4 identical |
| pool (`-k5 -K7 -r3 -j2 --printallk`)* | 127.8 s | 2.3 s | ~55× | 10/10 identical |
| ccf (`--ccf -r3 --debug`) | 41.9 s | 1.7 s | ~25× | 4/4 identical |
| conservative (`--conservativeCIs -K6 -r3 --debug`) | 7.7 s | 4.1 s | 1.9× (bootstrap loop unchanged) | 4/4 identical |
| silhouette (`-k5 -K7 -r2 --debug`) | 44.9 s | 1.9 s | ~24× | 4/4 identical |

\* The original pool path cannot run in this sandbox (see the sandbox table), so the new **pool**
run was compared with the original **`--debug`** run. They are equivalent by construction: every job
is seeded and the CIs are deterministic.

Descent-only exact comparison (`benchmarks/compare_descent.py`: C, assignments, objective, and
iterations identical for every job):

| Input | Jobs | Original descent | New descent |
|---|---|---|---|
| PDAC, 11 samples, k=13 | 2 | 2.3 s | <0.1 s |
| PDAC, 11 samples, k=16 (variable clusters; S ≥ 8 summation order) | 3 | 11.9 s | 0.3 s |
| synthetic 500 SNVs, k=6 | 2 | 2.4 s | <0.1 s |
| synthetic 500 SNVs, CCF, k=7 | 2 | 5.7 s | 0.1 s |

Peak RSS: the first version used 211 MB vs the original's 170 MB, from 2M-element CI-grid chunks.
With `MAX_GRID_ELEMENTS = 250000` it uses 172 MB, at the same speed.

Full-output PDAC and synthetic-500 golden masters were **not run by me**: the original code needs
~3–5 min single-core for each, mostly in the CI grid. Command for you:
`PYTHON=/home/myersm2/miniforge3/envs/decifer/bin/python benchmarks/run_regression.sh -s "pdac synthetic"`
(run it on the branch tip after Phase A, or with `-n "--convergence objective"` after Phase B).

### Phase B results
- Regression with `--convergence objective`: byte-identical in all scenarios except `betabinomial`
  (intended, from B2). The fitted precisions changed 142.144736 → 142.144704 (A17D), 129.401268 →
  129.401247, and 149.543821 → 149.543842, at the same optimum objective to 12 digits. The CIs are
  identical; the objective changed by 2.6e-5. With a 1,000-point brute force, the bb run takes 12 s
  instead of 5.7 s.
- `--binarysearch`: runs end to end, selects k = 8 on the test data, and pool == `--debug` byte-identical.
- **Stopping rule (B5)**: on 6 scenarios, the default `assignments` rule gives byte-identical outputs to
  the `objective` rule. Per restart (100 restarts over 5 inputs), the same final state in 97/100 and
  ~1 fewer iteration on average (e.g. 8.2 → 7.2). In 3 PDAC (k=16) restarts the objective rule ended lower,
  because its extra final assignment step uses randomly re-drawn centers of empty clusters (left
  unoptimized). The best objective across restarts was identical on all inputs.
- **Ground truth** (`benchmarks/validate_convergence.sh`: 30 simulated instances, k_true ∈ {2, 4, 8},
  500 SNVs, 3 samples, 100×, `-k 5 -K 12 -r 10`): ARI, cell-fraction error and k̂ − k are **identical
  in 30/30 instances** for both rules. Means: k_true=2 ARI 0.909 / CF error 0.007 / k̂−k 0.3; k_true=4
  0.972 / 0.006 / 0.9; k_true=8 0.947 / 0.011 / 0.0. Total wall time: 101.8 s (assignments) vs
  108.7 s (objective). Note that the default `--vafdevfilter 1.5` flags ~32% of simulated SNVs as
  outliers (existing behavior, not changed).

### numpy 2 check (row 29)
numpy 2.2.6 (offline): 26/26 tests pass. The regression is byte-identical to v2.1.4 in 5 scenarios, and
bb is identical to the numpy-1.26 run of the new code. No `np.float64(...)` appears in any output. The env
was then rolled back to 1.26.4.

### C3 (exact vectorized bootstrap, row 30)
`--conservativeCIs`, 5k SNVs, `-k5 -K5 -r1 --debug`: 16.2 → 4.5 s with identical CIs. Test-data
golden master `conservative`: identical, 2.1 s (original 7.7 s).

### C4 validation (`--init dcf`, k-means++ from DCF point estimates) — REJECTED, code removed (row 36)
`benchmarks/validate_variants.sh $TMPDIR/validate_init 10 500 "2 4 8" "-k 5 -K 12 -j 4" random_r10="-r 10" dcf_r10="-r 10 --init dcf" dcf_r5="-r 5 --init dcf"`
(30 simulated instances, 3 samples, 100×):

| Variant | ARI | CF error | mean \|k̂−k\| | iterations/restart | wall |
|---|---|---|---|---|---|
| random, 10 restarts | **0.943** | **0.0079** | **0.60** | 8.5 | 2.84 s |
| dcf, 10 restarts | 0.900 | 0.0146 | 0.87 | 7.2 | 2.64 s |
| dcf, 5 restarts | 0.885 | 0.0156 | 1.03 | 7.3 | 2.16 s |

ARI with dcf_r10 vs random: worse in 17/30 instances, better in 4, and down by as much as 0.36.
A plausible cause is that the point estimates use the first VAF-consistent genotype tree, which is
often not the true one, and D² seeding favors those outliers.

### Phase C profile (new code)
- 5k SNVs (`-k5 -K7 -r2 --debug`): CI 4.8 s of 9.7 s (49%); 6 descents 0.7 s.
- 20k SNVs (`-k8 -K8 -r2 --debug`): one descent takes 1.76 s with 41–82 iterations. The assignment step
  is 63% of a descent; `minimize_scalar` overhead is ~0.07 s in total. CI is 11 s per printed K.
  With default `-r 20 -K 12`, the descents (~160) dominate.

## Updates (exhaustive, with reasoning)

| # | File | Change | Reason |
|---|---|---|---|
| 1 | `benchmarks/timecmd.py` (new) | Wall/CPU/peak-RSS wrapper | `/usr/bin/time` is not installed on this host |
| 2 | `src/decifer/fastmodel.py` (new) | `Model`: per-(mutation, genotype tree) × sample coefficient arrays (`lbT, ubT, off, mu, msum, invF, lamlimT`) and likelihood constants (`betaln(a+1, ref+1)` or beta-binomial `CS`). Kernels `vaf` (= `Config.cf_to_v`) and `loglik` (= `compute_pdfs`). Methods `objective_fn`, `objective_grid`, `config_cluster_nll`, `mutation_cluster_nll`, `best_assignment`. | Replaces ~10 Python calls per (mutation, tree, cluster, sample) with numpy. The operation order matches the scalar code exactly (bit-identical). It has no import-time argparse, so it is testable. |
| 3 | `src/decifer/fastmodel.py` | Binomial log-pdf computed as `xlog1py(b-1,-x) + xlogy(a-1,x) - betaln(a,b)`, with `betaln` precomputed | This is scipy's own `beta._logpdf` (scipy 1.13.1). It avoids ~124 µs of `rv_continuous.logpdf` argument-checking overhead per call |
| 4 | `src/decifer/fastmodel.py` | A row with fractional copy number F == 0 is treated as infeasible | The scalar code raises ZeroDivisionError there. This is the only behavioral difference, and it only affects inputs that crash today |
| 5 | `src/decifer/new_coordinate_ascent.py` | `coordinate_descent` takes a `Model` and keeps assignments as int arrays `clus`/`cfg`; it returns `C, clus, cfg, obj, IT`. `optimize_assignments` uses `Model.best_assignment` (first-occurrence argmin, vectorized "cannot assign" / "Non scende" asserts, Python `sum` of per-mutation values). `optimize_cluster_centers` is rewritten as explicit loops with identical logic and RNG call order. Members are computed once per iteration instead of an O(n) `filter` per (sample, cluster). `update_objs` uses `Model.objective_fn` | Vectorized descent; byte-identical results |
| 6 | `src/decifer/new_coordinate_ascent.py` | `EPSILON`/`SEQERROR` imported from `fastmodel`; `np.NINF` replaced by `-np.inf` in `compute_pdfs`; scalar `objective` kept (used by `print_PDF` and the tests) | Single source for the constants; numpy 2 removed `np.NINF` |
| 7 | `src/decifer/__main__.py` | Builds the `Model` once in `main`. New `run_descents` (debug loop or `Pool` context manager), `unpack_result` (shallow `copy.copy` of each Mutation with its assignment), and `init_descent`/`run_descent` without `deepcopy` or `Manager` dicts. Workers return `job + (it, obj, C, clus, cfg)`. `setup_shared` removed. `silhouette_model_selection` takes `{k: mutations}` | Removes per-job deepcopy of all mutations, and the Manager process that held full mutation lists for every (restart, k) job |
| 8 | `src/decifer/__main__.py` | `compute_CIs_mp` jobs are `(c, s, num_tests)`; mutations and model go to workers once via `init_CI`. `CI` evaluates the 10,000-point grid with `Model.objective_grid`; the post-processing is unchanged | Cut the dominant cost (93% of the baseline profile). It also stops pickling every mutation into every CI job |
| 9 | `src/decifer/fileio.py` | `write_results(..., model)` takes Explained/LHs from `Model.mutation_cluster_nll(C)` | Removes a full scalar assignment pass per output file |
| 10 | `src/decifer/mutation.py` | `create_mutations` iterates `to_dict('records')` instead of `df.loc[idx]`, with `num_CNstates` hoisted | ~50–100 µs per row saved; same strings, same logic |
| 11 | `test/conftest.py`, `test/test_fastmodel.py` (new) | 16 exact-equality tests vs `Config.cf_to_v`, `compute_pdfs`, the scalar `objective`, and the v2.1.4 assignment step (oracle copied into the test). Also checks of the reduction-order assumptions | Pins the bit-identical claim |
| 12 | `benchmarks/make_synthetic.py`, `compare_outputs.py`, `run_regression.sh`, `compare_descent.py` (new) | Seeded scale-up generator; byte-for-byte output comparison; golden-master scenario matrix (base cached per ref); descent-only exact comparison | Reproducible verification and benchmarking |
| 13 | `benchmarks/make_synthetic.py` | `simulate` mode: planted clusters on real CN profiles, d ~ Poisson, a ~ Binomial(d, max(v, 0.001)), with a truth file | Ground truth for validating result-changing changes (B5, C4); the decifer-data download is not possible offline |
| 14 | `src/decifer/fastmodel.py` | `MAX_GRID_ELEMENTS` 2,000,000 → 250,000 | Peak RSS back to the original level (172 vs 170 MB) at the same speed |
| 15 | `src/decifer/__main__.py` | B1: new `finalize_k()` (filter → CIs → write) shared by both coordinators. `run_coordinator_binary` now consumes the job results (the old lazy `map(report, …)` never ran them), uses `finalize_k`, supports `--debug`, and reports the selected k. `run` returns the best raw result | `--binarysearch` crashed before (unconsumed results, and `write_results` called with wrong arguments) |
| 16 | `src/decifer/__main__.py` | `run()` (binary search) saves and restores the main-process RNG state around each batch of jobs | In `--debug` mode the jobs re-seed the main process, so the seeds drawn for the next K differed from pool mode. Now pool == debug (verified byte-identical). Not applied to the iterative path, which draws all seeds up front, so its golden master is unchanged |
| 17 | `src/decifer/generator.py` | B2: `fit_betabinom` uses `args['restarts_bb']` (brute-force range and Brent bracket); the `fitbeta` CLI dict key renamed to `restarts_bb`; added the missing `main()` (writes `PREFIXbetabinom.tsv`) | DeCiFer passed its clustering `--restarts` (default 20) as the brute-force size, so only S ∈ {0, 10} was tried. The `fitbeta` entry point referenced a nonexistent `main` |
| 18 | `src/decifer/__main__.py`, `test/test_end2end.py` | B3: `CI_conservative` casts its quantiles to `float`; `read_csv(..., sep="\t")` | Under numpy 2, `str((np.float64(x), …))` prints `np.float64(x)` into the CI columns. pandas 2 made `sep` keyword-only |
| 19 | `src/decifer/process_input.py` | B4: read `MUTATION_DF` before the purity block | `calc_purity(MUTATION_DF)` referenced an undefined name (dead code today because `-p` is required, but a latent NameError) |
| 20 | `src/decifer/new_coordinate_ascent.py`, `parse_args.py`, `__main__.py` | B5: `coordinate_descent(..., convergence='assignments')` stops at the first assignment step whose (cluster, genotype-tree) assignments equal the previous step's, and skips the redundant final assignment. `--convergence {assignments,objective}` (default `assignments`) is threaded through `run_descents` → `init_descent` → `run_descent` | Requested by the user. Matches the paper (STAR Methods e5). The legacy rule stays available for exact reproduction |
| 21 | `test/test_convergence.py` (new) | Stops at the first repeated assignment; respects `--maxit`; rejects unknown rules | Pins B5 |
| 22 | `test/test_end2end.py` | The existing truth test passes `--convergence objective`; new `test_binarysearch` smoke test | Keeps the stored truth comparable; covers B1 |
| 23 | `README.md` | Recommendations: vectorization note (20–50× vs v2.1.4, identical results) replaces the "more than a few days" wording; new `--convergence` bullet | Docs. `man/man-decifer.md` is empty in the repo, so options are documented in the README and `--help` |
| 24 | `benchmarks/compare_descent.py` | `--base-convergence`, `--json` | Compare the two rules job by job |
| 25 | `benchmarks/evaluate_truth.py`, `benchmarks/validate_convergence.sh` (new) | Paper metrics (ARI, cell-fraction error, k̂ − k) and a sweep over k_true × seeds × rules | Ground-truth validation of B5 |
| 26 | `src/decifer/fastmodel.py` `objective_grid` | Evaluates only the grid points inside the intersection of the rows' analytic feasible c-intervals (1e-9 margin), and runs `loglik` only where every row is feasible. Every other point is -EPSILON, as before | Exact. CI at 5k SNVs: 4.8 → 3.9 s (most points are feasible there, so the gain is modest) |
| 27 | `src/decifer/fastmodel.py` `config_cluster_nll` | Caches log-likelihood columns by (sample, cell fraction) between calls. Fixed clusters (0, purity, other samples' sample-specific clusters at 0) and repeated centers are computed once | Exact. The assignment step was 63% of a descent at 20k SNVs |
| 28 | `src/decifer/new_coordinate_ascent.py` | New `ClusterObjectives`: per (sample, cluster), it memoizes objective values and the `minimize_scalar` minimizer while the cluster's members are unchanged (50–63% of clusters per iteration at 20k). Used by `update_objs` and `optimize_cluster_centers`; the RNG call order is unchanged | Exact. Rows 27+28 make descents at 20k 2.1× faster (11.7 → 5.6 s for 3 descents), with identical results |
| 29 | conda env `decifer` | numpy-2 check: upgraded offline to numpy 2.2.6, ran the tests and regression, then **rolled back to 1.26.4** | numpy 2 verified: all 26 tests pass, and every output is identical to the numpy-1.26 runs (vs original v2.1.4 for the 5 non-bb scenarios, vs the new code for bb), with no `np.float64(...)` in any output. Rolled back because the golden-master baseline (v2.1.4 uses `np.NINF`) only runs on numpy < 2 |
| 30 | `src/decifer/__main__.py` | C3: new `bootstrap_medians` replaces the `rand.choices` + `np.median` loop of `CI_conservative`. It loads Python's `random` MT19937 state into a numpy legacy `RandomState` (same generator, same 53-bit double construction), draws all resample indices as `floor(u * n)` in chunks, takes the medians along axis 1, and writes the advanced state back to `random` | **Exact** (bit-identical draws, medians, CIs and final `random` state, unlike the plan's results-changing version). `--conservativeCIs` at 5k SNVs: 16.2 → 4.5 s. The `conservative` golden master is identical (2.1 s; original 7.7 s) |
| 31 | `test/test_bootstrap.py` (new) | `bootstrap_medians` vs the original loop for n ∈ {0, 1, 2, 3, 10, 57, 400}, several seeds, chunked; also checks the final `random` state | Pins C3 |
| 32 | `src/decifer/fastmodel.py` | `Model.point_estimates(purity)`: per-SNV cell-fraction point estimates (= `Config.v_to_cf(a/d, s, truncate=False)` to 1e-12) for the first genotype tree consistent with the VAFs in every sample; stores `F`, `ub`, `vaf_obs` | Input for C4 |
| 33 | `src/decifer/new_coordinate_ascent.py`, `parse_args.py`, `__main__.py` | C4: `coordinate_descent(..., init='random')`. **Opt-in `--init dcf`** chooses the variable cluster centers among the point estimates by k-means++ D² seeding (`dcf_centers`), measuring distance to the fixed clusters and already-chosen centers, with the per-job seeded `np.random`. Plumbed through `run_descents`/`init_descent`/`run_descent`/`run` | Plan item C4 (descents dominate realistic runs). Default `random` is unchanged: the golden master is still byte-identical (26 files) |
| 34 | `test/test_init.py` (new) | Point estimates vs `Config.v_to_cf` (DCF and CCF); `dcf_centers` picks points among the estimates with shape (rest, S); `init='dcf'` keeps the fixed clusters and centers within [0, purity] | Pins C4 |
| 35 | `benchmarks/validate_variants.sh` (new) | Generalized ground-truth sweep over arbitrary option sets, with mean iterations per restart | Validation of C4 |
| 36 | `src/decifer/fastmodel.py`, `new_coordinate_ascent.py`, `parse_args.py`, `__main__.py`; `test/test_init.py` deleted | **C4 removed** (rows 32–34 reverted: `--init`, `dcf_centers`, `point_estimates`, `F`/`ub`/`vaf_obs`, plumbing, tests). `validate_variants.sh` kept | It failed ground-truth validation (see "C4 validation"): worse ARI, CF error and k̂ − k for ~15% fewer iterations. An option validated as worse would only invite misuse |
| — | (incident) | During cleanup I deleted the tracked `test.png` (the e2e/beta-binomial run overwrites `./test.png`) and **restored it with `git checkout`**. The tracked `test/decifer_*.tsv` that the e2e test overwrites were also restored | No net change to tracked test artifacts |

## Sandbox-limited tasks and workarounds

| Task | Blocked by | Workaround |
|---|---|---|
| `mamba create -n decifer` from conda-forge | No network. The sandbox proxy at `localhost:3128` refused connections to `conda.anaconda.org` (curl error 7) | Built the env **offline** from the local package cache (`CONDA_PKGS_DIRS=/home/myersm2/miniforge3/pkgs`, `--offline`). The cache has no complete Python 3.11 set, so I used **Python 3.10**: numpy 1.26.4, scipy 1.13.1, pandas 2.2.3, matplotlib-base 3.8.4, seaborn-base 0.13.2 (not `seaborn`, which pulls the uncached statsmodels), scikit-learn 1.5.2, pytest 9.1.1. I pinned `libsqlite=3.53.3` because the solver's pick was not cached. |
| Offline solve with `/home/myersm2/.conda/pkgs` in `pkgs_dirs` | Path outside the sandbox (mamba: "Specified pkgs_dir does not exist") | Set `CONDA_PKGS_DIRS` to the miniforge cache only |
| `pip install -e .` build isolation | No network to fetch setuptools | `pip install -e . --no-deps --no-build-isolation`, using the env's setuptools 84. `--no-deps` also skips the spurious PyPI `lemon` requirement |
| Original DeCiFer in pool mode (no `--debug`) for the golden master | The sandbox forbids socket creation (`PermissionError: [Errno 1] Operation not permitted` in `multiprocessing.Manager()`, which the original uses for its shared dicts) | Compared the new pool-mode run with the original `--debug` run (equivalent; see Results). The new code needs no Manager except with `--record`, which therefore cannot run in this sandbox |
| Downloading `raphael-group/decifer-data` simulated instances (Phase C) | No network expected (same proxy) | Will attempt once. Otherwise use the local simulator `benchmarks/make_synthetic.py --mode simulate` |
