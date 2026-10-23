# 2026-10-21 — Faster, exact state-tree enumeration (`generatestatetrees`)

Approved plan: `/home/myersm2/.claude/plans/can-the-enumeration-of-frolicking-widget.md`
Branch: `main` (working tree; nothing committed. An earlier version of this log said
`perf/vectorize-core`, which was wrong: `git branch --show-current` gives `main`)
User decisions: target CN ≥ 8; output byte-identical to the current enumerator; C++ only.

## Summary
**Status (2026-10-22): the three requested items are done, verified and documented.** Nothing is
committed.

1. **Symmetry fix** (row 24) compiled and verified. Shipped file with `-symmetry -threads 3`:
   0.33–0.42 s (was ~110 s), byte-identical for tree, flashlight and auto.
2. **Dominator ("mandatory cell") relaxation** in the witness oracle (row 27). Exact (sound by
   induction, see row 27); byte-identical on every gate list. Clustered sets solved within 10 s:
   36/90 → 48/90 (|L| = 3: 16/30 → 28/30).
3. **Frontier DP** (rows 28, 30–38): a new exact algorithm, `-algorithm frontier`, also used
   - inside flashlight, as the oracle's exact fallback: an undecided prefix is solved by
     enumerating all its state trees with the constrained DP, so no oracle call is unbounded
     unless the DP exceeds its memory limit;
   - in `auto`, which now races tree, flashlight and frontier and is the **default**.

   Byte-identical to the v2.1.4 reference and to the other engines on every list (shipped, exh3,
   exh4, ladder, random, combinations mode). The constrained mode was checked with a build that
   sends every unrouted oracle call to the DP (~13,000 DP calls, all identical).

Within 10 s per set, the 90 clustered sets: 36 (start of the session) → 48 (dominator) → 54 (DP
fallback) → **57 (`auto`)**. Of the 33 unfinished, 31 have 4–5 states.

**Limits.** The DP is exponential in the lattice width, because paths can wiggle across the
frontier, and labels multiply its states. The longest CN 6 set I traced
({1,1;3,5;3,6;5,5}) now takes 22 s. A CN 8 label-phase DP ({1,1;6,7;6,8;7,8}) keeps ~640 k
states per step and took more than 300 s. Sets with |L| ≥ 4 at CN ≥ 8 with both alleles high
remain slow or out of reach. `-frontierLimit` (default 2 M states per step, about 1 GB) bounds
the DP's memory; over the limit, flashlight falls back to its exhaustive search.

Tests: `test/test_statetrees.py` 27 passed. Docs: `scripts/README.md`. Commands for longer runs
(for the user): see "Commands for the user" below.

## Method
- Exactness target: Out(L) = {S(T) : T truly valid} (see plan). Every pruning rule must be
  monotone, so it only removes branches that contain no valid T. The output goes through the
  unchanged `StateGraph::_dict` / `operator<<`, so equal sets give byte-identical files.
- Reference: `generatestatetrees` built from `git archive 1b6f6d8 src/decifer/cpp`, in
  `$TMPDIR/st-ref/build` (`$TMPDIR` = `/tmp/claude-164077095`).
- Build (both binaries): system `/usr/bin/cmake` 3.26.5 and `/usr/bin/g++` 8.5,
  `-DLIBLEMON_ROOT=/home/myersm2/miniforge3/pkgs/lemon-1.3.1-h9457a4f_11` (headers only),
  `-DBOOST_ROOT=/home/myersm2/miniforge3/envs/hatchet` (Boost 1.75, read-only use). The new
  build is in `$TMPDIR/st-new` (built with `-DCMAKE_CXX_FLAGS=-Wall`).
- Set lists in `$TMPDIR/st-sets/` (from `benchmarks/statetrees/make_cn_sets.py`):
  `shipped.txt` (469), `exh4.txt` (2,325 = every L ∋ (1,1), |L| ≤ 4, CN ≤ 4), `exh3.txt` (576),
  `ladder.txt` (45: families diag, diag3, thin1, thin2, wgd for k = 2..10).
- Runs in `$TMPDIR/st-runs/<name>/` (from `benchmarks/statetrees/run_sets.py`; `runs.jsonl` per run).
- Build since 2026-10-22: `$TMPDIR/st-dev`, `-DLIBLEMON_ROOT=/home/myersm2/miniforge3/envs/decifer-cpp`
  (that env's Boost 1.85 headers now come first on the include path).

### Frontier DP (`-algorithm frontier`)
- A truly valid tree T is, on the lattice, a choice of used cells, the label of each used cell's
  node (the mutation cell holds the unmutated node, label 0, and the mutated node, label 1) and
  arcs between adjacent used cells. Cells that cannot reach L are dropped.
- Cells are processed in rows along the shorter side of the lattice (frontier width
  min(X, Y) + 1). Per frontier cell the state keeps: used, mutation cell, has parent, has child,
  label, owner (index in L, *unknown* = above its component's top, or *none* = below a root not in
  L) and component. Per component: the nodes of L whose state tree parent is the owner of the
  component's future parent (*pending*). Plus: mutation placed, component of the root.
- Decisions per cell: unused (not for L or the root), a node with label 0..x (maternal run) or
  0..y (paternal run), or the mutation (cell in L, once). Per frontier neighbor: no arc, or an arc
  in either direction, if the state graph has it for the two labels (mutation cell: from either
  node, into the unmutated node), into a node without a parent, between different components
  (so T stays a forest). Attaching a component resolves its *unknown* owners and pending nodes.
- A cell leaving the frontier must have a parent (root excepted) and, if not in L, a child. When
  the root's component leaves the frontier, T is complete if nothing else is on the frontier,
  every cell of L is processed and the mutation is placed.
- Each state keeps the set of partial state trees (parent and label of every processed node of L)
  of all partial T reaching it: the transitions depend only on the state, so equal states have
  equal completions. Accepted partial state trees of both allele runs are the output.
- Exactness: accepted decision sequences are exactly the truly valid trees (one sequence each), and
  the state tree is read off the owners, as in `finalize`.

## Progress
- [x] Phase 0: build, reference build, set generators, runner, golden master on the shipped sets
  - reference on shipped (5 s timeout, 4 jobs): 304/469 finished, 165 timed out
  - [x] reference on `exh3.txt` (10 s timeout, 2 jobs) finished after the pause: 488/576 sets
    finished, 88 timed out, 990 s total wall (`$TMPDIR/st-runs/ref-exh3`). **Not yet compared**:
    on resume, run tree on `exh3.txt` and compare (commands below)
- [x] Phase 1: `-algorithm tree` (pruned engine) — gate passed (see Results)
- [x] after resume: tree on `exh3.txt` vs reference: 488/488 identical (the 88 others timed out in
  the reference); 990 s (reference) vs 4.4 s (tree)
- [x] Phase 2: cache check, `-threads`, `-symmetry` (fixed and verified 2026-10-22)
- [x] Phase 3: `-algorithm flashlight` (gate passed), dominator relaxation, DP fallback
- [x] Frontier DP: `-algorithm frontier`, constrained mode in flashlight, three-way `auto` (default)
- [x] Phase 4: tests (27 pass), `scripts/README.md`, benchmark set lists in
  `benchmarks/statetrees/sets/`
- [ ] For the user (expensive): longer-timeout benchmark at CN 8–10, and DP vs flashlight
  cross-validation where both finish (commands below)

### To resume (from the pause; done)
```bash
cd /data1/shahs3/users/myersm2/repos/decifer
PY=/home/myersm2/miniforge3/envs/decifer/bin/python
T=/tmp/claude-164077095
# rebuild the new engine (if $TMPDIR was cleaned, re-run cmake as under Method first)
make -C $T/st-new generatestatetrees
# check the reference exhaustive run that was still going when paused, against a tree run on the same list
$PY benchmarks/statetrees/run_sets.py run $T/st-new/generatestatetrees $T/st-sets/exh3.txt $T/st-runs/tree-exh3 --jobs 2 --args=-stats
$PY benchmarks/statetrees/run_sets.py compare $T/st-runs/ref-exh3 $T/st-runs/tree-exh3   # expect 0 DIFFERENT
```

## Results

### Reference (v2.1.4) timings
- {(1,1),(3,3)}: 5.6 s; {(1,1),(4,4)}: > 120 s (timeout); {(1,1),(2,2),(2,3),(3,2)}: 5.5 s
  (single runs, loaded node).
- Shipped sets, 5 s timeout, 4 parallel jobs: 304 finished (39.3 s total), 165 timed out.

### Phase 1 gate — PASSED
- `generatestatetrees src-sets/shipped.txt` (all 469 sets in one process): output is
  **byte-identical to the shipped `src/decifer/state_trees.txt`**, in 1.3 s wall.
- Per set vs the reference: 304/304 identical (the sets the reference finished); total wall
  39.3 s (reference) vs 2.4 s (tree, mostly process start-up).
- Small spot checks, identical: {0,0;1,0;1,1} (12 trees, 49 search calls), {1,1}, {1,1;2,2},
  {1,1;3,3} (423 calls vs > 10⁶ in the reference), {1,1;2,2;2,3;3,2} (228 trees, 3,395 calls),
  {0,0;1,1}.
- `exh4.txt` (2,325 sets, CN ≤ 4) with tree: all finished, 581 s total wall with 2 jobs
  (not yet compared: no reference run on exh4).

### Ladder (tree engine, 60 s timeout, 2 jobs)

| family | k: wall (trees) |
|---|---|
| diag {(1,1),(k,k)} | 2–4: ≤ 6 ms (3); 5: 0.15 s (3); 6: 20.8 s (3); 7–10: timeout |
| diag3 {(1,1),(2,2),(k,k−1)} | 2: 40 ms (22); 4: 0.21 s (26); 5: 0.79 s (27); 6: 56 s (28); 7–10: timeout |
| thin1 {(1,1),(k,1)} | 2–10: ≤ 8 ms (4) |
| thin2 {(1,1),(k,2),(k,0)} | 2–10: ≤ 0.7 s (15–24) |
| wgd {(1,1),(2,2),(2,k),(k,2)} | 3: 16 ms (228); 4: 0.15 s (306); 5: 35 s (364); 6–10: timeout |

Note: the number of state trees grows slowly with k (diag: always 3), while the search grows with
the number of self-avoiding paths. This motivates the S-level (flashlight) algorithm.

### Phase 3 gate (flashlight vs tree; all byte-identical)
- shipped file (469 sets): reproduced byte for byte (now 0.27 s with `-threads 2`)
- `exh3.txt` 576/576; `exh4.txt` 2,325/2,325 (run with an earlier flashlight build)
- all-combinations mode `-maxXY 3 -maxCN 3` (560 sets) and `-maxCN 4` (1,820 sets, including sets
  without (1,1)); `-maxXY 2 -maxCN 2` also identical to the reference and legacy
- ladder: all 32 sets both finished
- `-threads 3/4` identical to sequential for tree and flashlight
- `pytest test/test_statetrees.py`: 10 passed

### Phase 3 tuning (oracle = witness search), timings on a loaded node

| set | tree | flashlight v1 | + parent-not-in-T check | + label order, arc heuristic | + candidate propagation | + greedy routing |
|---|---|---|---|---|---|---|
| {0,0;0,4;1,1;4,0} | 8.5 s | 34.6 s | 0.12 s | 0.013 s | 0.0007 s | |
| {0,0;0,4;1,1;4,4} | 6.5 s | | 1.6 s | 0.24 s | 0.057 s | 0.050 s |
| {1,1;2,2;6,5} | 56 s | | | 0.67 s | 0.0012 s | |
| {1,1;2,2;7,6} | > 60 s | | | > 120 s | 0.0022 s | 0.0008 s |
| {1,1;2,2;2,5;5,2} | 35 s | | | 3.3 s | 1.57 s | 0.38 s |
| {1,1;2,2;2,6;6,2} | > 60 s | | | > 120 s | > 120 s | 70 s |
| {1,1;(k,k)}, k ≤ 10 | k=6: 20.8 s, k ≥ 7: > 60 s | ≤ 16 ms | | | | |

What each step fixed (found by instrumenting a scratch copy of the code in `$TMPDIR/st-dbg`
that prints every oracle call):
1. parent-not-in-T check: prefixes such as "parent of (0,4) is (0,0)" ((0,0) is a sink) were
   only refuted after millions of partial trees, because the check needed the parent in T.
2. labels decided parents-first, and the witness search picks the frontier arc closest to an
   allowed target (constraint- and label-aware distances over genotype nodes).
3. candidate propagation: e.g. "(7,6) label 6 below a mutation in (1,1), (2,2) label 0 below
   (7,6)" is infeasible (label 6 at x=7 needs the paternal allele, whose copies then all stay
   mutated). The check started from every node (7,6) could be, including unreachable ones.
4. greedy routing: most remaining time was in feasible prefixes whose witness the exhaustive
   search found only after a long search.

5. routing strategies (shortest-first, longest-first, 8 random orders with random tie-breaking):
   {1,1;2,5;10,2} timeout → 8 ms (greedy shortest-first walled off the only route).
6. negotiated congestion routing (PathFinder) when all parents are decided, and a root-state
   relaxation check before any routing (infeasible prefixes return at once).

**Limitation found** (step 6 analysis): some prefixes are infeasible for *topological* reasons
that no per-edge relaxation sees. Example: L = {1,1; 1,10; 2,0; 7,2} (max x = 7), with (1,10) and
(7,2) children of (1,1) and (2,0) a child of (1,10). The paths (1,1)→(7,2) and (1,1)→(1,10)
enclose rows y ≤ 1, so (2,0) cannot be reached from (1,10) (a planar crossing obstruction). The
exhaustive search refutes this only after trying every shape of the (1,1)→(7,2) path, which is
exponential. An exact method that handles this needs topological reasoning, e.g. a
frontier-based (transfer-matrix) dynamic program over the lattice.

### Random sets, frozen binary `gst-snapshot1` (30 s timeout, 3 jobs; finished/10 per cell)

| set kind | CN | \|L\|=3 flashlight / tree | \|L\|=4 | \|L\|=5 |
|---|---|---|---|---|
| random "thin" (each state has one allele ≤ 2) | 6 | 10 / 9 | 9 / 7 | 5 / 4 |
| | 8 | 10 / 7 | 4 / 3 | 2 / 1 |
| | 10 | 10 / 4 | 4 / 3 | 0 / 0 |
| random square | 6 | 10 / 8 | 6 / 5 | 1 / 1 |
| | 8 | 8 / 2 | 2 / 2 | 0 / 0 |
| | 10 | 9 / 0 | 0 / 0 | 0 / 0 |
| clustered (flashlight only) | 6 | 8 | 6 | 4 |
| | 8 | 5 | 5 | 3 |
| | 10 | 3 | 2 | 0 |

Every set tree finished, flashlight finished too, with identical output (38/38 thin, 18/18 square).
Example failure analysis, clustered {1,1;6,5;6,6}: (6,6)'s only entry is (5,6) ((6,5) is in L and
the lattice ends at 6), so a prefix where (6,5) needs label 2 below the mutated copy of (6,6) is
infeasible (that route also needs (5,6)). The relaxations don't reserve (5,6) for the path into
(6,6), so the search tries every such path before refuting.

### 2026-10-22: symmetry fix, dominator relaxation, frontier DP

**Symmetry fix** (shipped file, all 469 sets in one run; byte-identical to `state_trees.txt` in
every case):

| algorithm | 1 thread | 3 threads | `-symmetry` | `-symmetry -threads 3` |
|---|---|---|---|---|
| tree | 1.31 s | 0.51 s | 0.81 s | 0.34 s (was ~110 s) |
| flashlight | 1.77 s | 0.66 s | 1.05 s | 0.42 s |
| auto (tree + flashlight) | 0.82 s | 0.49 s | 0.52 s | 0.33 s |

**Gates of the new code** (`run_sets.py check`, one process per list; every set byte-identical):

| list (sets) | flashlight + dominator | frontier (enumeration) | flashlight with `GROW_BUDGET = 1` (DP for every unrouted call) |
|---|---|---|---|
| shipped (469) | ✓ | ✓ (and reproduces `state_trees.txt`) | ✓ (4,011 DP calls; with `-frontierLimit 50`: ✓) |
| exh3 (576) | ✓ | ✓ | ✓ (4,946 DP calls; with `-frontierLimit 50`: ✓, 2,366 over the limit) |
| exh4 (2,325) | ✓ (13 s; tree 581 s) | ✓ (422 s, before rows 35–36) | ✓ (337 s) |
| ladder (32 tree-finished) | ✓ | ✓ | |
| random square (18) / thin (38) | ✓ / ✓ | ✓ (15 s / 1.0 s after row 36; tree 58.9 s / 21.9 s) | |
| combinations `-maxXY 3 -maxCN 3` / `-maxCN 4` | ✓ / ✓ | ✓ / ✓ | |

`auto` (tree + flashlight + frontier) on random thin / square: 0.5 s / 2.1 s.

**Frontier DP scaling** (enumeration, both alleles; after rows 35–36):

| set | lattice used | states (max per step) | time |
|---|---|---|---|
| {1,1;4,4} | 4×4 | 69 k (7.4 k) | 0.09 s |
| {1,1;5,5} | 5×5 | 1.6 M (158 k) | 1.9 s |
| {0,0;0,4;1,1;4,4} | 5×5 | 1.9 M (252 k) | 2.5 s |
| {1,1;2,2;5,4} | 5×5 | 0.44 M (65 k) | 0.5 s |
| {1,1;6,6} | 6×6 | | > 60 s (before rows 35–36) |

Labels multiply the states about 10–30× (instrumented: e.g. 99,670 states vs 4,820 without
labels), mostly in path fragments not yet attached to the root. In constrained mode the path
"dips" (segments that go down and come back) dominate: a CN-8 call ({1,1;6,7;6,8;7,8}, all
parents decided, 7-wide frontier) holds ~640 k states per step.

**Flashlight with the DP fallback**, {1,1;3,5;3,6;5,5} (CN 6, |L| = 4; never finished before):
104 s (DP as existence oracle) → 117 s (DP enumeration per undecided decision) → 42 s (one DP for
all labels per parent configuration) → 25 s (row 35) → 22 s (row 36); 471 trees, identical output
in every version.

**Clustered sets** (`benchmarks/statetrees/sets/clustered.txt`), sets finished within the timeout
(per cell out of 10); every pair of runs that both finished a set gave identical output:

| clustered | flashlight v1 (30 s) | + dominator (10 s) | + DP fallback (10 s) | after rows 35–37 (10 s) | `auto` (10 s) |
|---|---|---|---|---|---|
| CN 6, \|L\| = 3 / 4 / 5 | 8 / 6 / 4 | 10 / 6 / 4 | 10 / 6 / 6 | 10 / 6 / 6 | 10 / 6 / 6 |
| CN 8, \|L\| = 3 / 4 / 5 | 5 / 5 / 3 | 10 / 5 / 3 | 10 / 6 / 3 | 10 / 6 / 3 | 10 / 6 / 4 |
| CN 10, \|L\| = 3 / 4 / 5 | 3 / 2 / 0 | 8 / 2 / 0 | 8 / 3 / 2 | 8 / 3 / 2 | 8 / 4 / 3 |
| total | 36 | 48 | 54 | 54 | **57** |

On the 54 sets both finished, `auto` took 8.9 s in total, flashlight 34.8 s.

**Final gate** (final build `$TMPDIR/st-dev/gst-final`, default algorithm `auto`, `-threads 2`):
shipped file reproduced byte for byte; tree-exh3 576/576, tree-exh4 2,325/2,325 (33 s; tree 581 s),
tree-ladder 32/32, random square 18/18, random thin 38/38, v2.1.4 reference on exh3 488/488 (0.5 s;
reference 990 s) and shipped 304/304; combinations mode `-maxXY 3 -maxCN 3` and `-maxCN 4`
identical; pytest 27 passed. A clean configure against `decifer-cpp` only
(`-DLIBLEMON_ROOT=$E -DBOOST_ROOT=$E`) builds and reproduces the shipped file.

## Commands for the user

Expensive runs I did not make. Build into the repository's `build/` (the `$TMPDIR` of my sandbox is
node-local):

```bash
cd /data1/shahs3/users/myersm2/repos/decifer
E=/home/myersm2/miniforge3/envs/decifer-cpp
PY=/home/myersm2/miniforge3/envs/decifer/bin/python
mkdir -p build && (cd build && cmake ../src/decifer/cpp -DLIBLEMON_ROOT=$E -DBOOST_ROOT=$E -DCMAKE_BUILD_TYPE=Release && make generatestatetrees)
GENERATESTATETREES=build/generatestatetrees $PY -m pytest -q test/test_statetrees.py   # expect 27 passed

# 1. how many hard sets finish within 10 minutes (auto = 3 threads per set; 4 sets at a time)
O=benchmarks/statetrees/out; mkdir -p $O
sbatch -c 12 --mem 16G -t 12:00:00 -o $O/auto-600.slurm.log --wrap \
  "$PY benchmarks/statetrees/run_sets.py run build/generatestatetrees benchmarks/statetrees/sets/clustered.txt $O/auto-600-clustered --timeout 600 --jobs 4 --args='-algorithm auto -stats'"

# 2. cross-validation at CN >= 8 of the two independent exact methods (frontier DP without a
#    state limit vs flashlight), on the sets both finish
sbatch -c 4 --mem 32G -t 12:00:00 -o $O/xval.slurm.log --wrap \
  "$PY benchmarks/statetrees/run_sets.py run build/generatestatetrees benchmarks/statetrees/sets/clustered.txt $O/frontier-clustered --timeout 1800 --jobs 2 --args='-algorithm frontier -frontierLimit 0 -stats' && \
   $PY benchmarks/statetrees/run_sets.py run build/generatestatetrees benchmarks/statetrees/sets/clustered.txt $O/flashlight-clustered --timeout 1800 --jobs 2 --args='-algorithm flashlight -stats' && \
   $PY benchmarks/statetrees/run_sets.py compare $O/frontier-clustered $O/flashlight-clustered"
```

Verification and follow-up:
- Run 1: `tail -1 $O/auto-600.slurm.log` prints `90 sets, N finished`; N should be ≥ 57.
  `grep -h '"exit": 124' $O/auto-600-clustered/runs.jsonl` lists the sets that still do not finish.
- Run 2: the last line must say `0 DIFFERENT`. A difference would mean a bug in one of the two
  methods (they share no search code). The `frontier` runs use no state limit, so watch memory
  (`sacct -j <id> -o MaxRSS`); lower `--jobs` if a job is killed for memory.
- If Run 1 leaves important sets unfinished, the next step is reducing the label dimension of the
  DP (see Results: labels multiply its states 10–30×).

## Updates (exhaustive, with reasoning)

| # | File | Change | Reason |
|---|---|---|---|
| 1 | `ai-logs/2026-10-21-statetree-enumeration.md` (new) | This log | User instructions |
| 2 | `benchmarks/statetrees/make_cn_sets.py` (new) | Writes set lists: `shipped` (from a state-tree file), `exhaustive`, `random` (square/thin, seeded), `ladder` (diag, diag3, thin1, thin2, wgd) | Reproducible benchmark and golden-master inputs |
| 3 | `benchmarks/statetrees/run_sets.py` (new) | `run`: one process per set with a timeout, per-set outputs and `runs.jsonl`, resumable, `--jobs`. `compare`: byte comparison of two runs on the sets both finished | Per-set timeouts are needed because the reference does not finish many sets. Written in Python instead of the planned `run_sets.sh` (simpler timeouts and JSON) |
| 4 | `src/decifer/cpp/enumerator.{h,cpp}` (new) | `StateTreeEnumerator`: flat-array reimplementation of the Gabow–Myers frontier search of `StateGraph::grow`, with rules R1 (mutation arcs only from L states), R2 (drop nodes whose copy state cannot reach L), R3 (true leaf bookkeeping; stop once L is covered and the mutation placed), R4 (non-L leaves ≤ uncovered L states), R5 (lattice BFS reachability), R6 (mutation still placeable). S is built from incrementally maintained owners (nearest L ancestor) | Phase 1. Exact (monotone pruning; see plan) and verified byte-identical |
| 5 | `src/decifer/cpp/enumerator.{h,cpp}` | Witness mode: `Constraints` (mutation state, parents, labels), `findWitness`, `consistent` (checked when an L node is attached), `reachableWitness`/`reachesChild` (node-level BFS through unused non-L copy states, honoring the parent/label constraints), `_nodeOfIndex` | Oracle for Phase 3. **Not yet compiled/tested** |
| 6 | `src/decifer/cpp/stategraph.{h,cpp}` | `hasStateTrees`, `addStateTrees`; `getStateTrees` checks `_dict` before building the state graph | Dispatch from main without changing the output path; avoid building G for cached sets |
| 7 | `src/decifer/cpp/generatestatetreesmain.cpp` | `computeStateTrees` dispatch; options `-algorithm tree\|legacy` (default `tree`) and `-stats` (per-set search statistics and time on stderr) | Select the engine; keep the original for comparison |
| 8 | `src/decifer/cpp/CMakeLists.txt` | `enumerator.cpp`/`.h` added to `generatestatetrees` | Build |
| 9 | `src/decifer/cpp/enumerator.cpp` | Removed a vacuous check in `consistent` (mutation state attached after the mutation) | It could never fire |
| 10 | `src/decifer/cpp/flashlight.{h,cpp}` (new) | `FlashlightEnumerator`: decisions mutation state → parents (index order; cycle check; -2 "no parent" only if (1,1) ∉ L) → labels of descendants of the mutated copy; oracle = `StateTreeEnumerator::findWitness`; a decision equal to the current witness needs no oracle call; at the leaves the witness's S is inserted | Phase 3: enumerate state trees without enumerating paths |
| 11 | `generatestatetreesmain.cpp`, `CMakeLists.txt` | `-algorithm flashlight`; flashlight statistics in `-stats`; flashlight sources in CMake | Phase 3 |
| 12 | `generatestatetreesmain.cpp` | `-threads N`: `computeAllStateTrees` collects the sets of each input file (or of the all-combinations mode) and computes them with a thread pool (`enumerateStateTrees` is thread safe); results go to the ordered `_dict` in input order; with 1 thread (or `legacy`) the original sequential progress output is kept | Phase 2 |
| 13 | `test/test_statetrees.py` (new) | Shipped-file reproduction (tree/flashlight × 1/3 threads), agreement with legacy (CN ≤ 2 sets, all-combinations mode), known trees, tree vs flashlight up to CN 5, flashlight at CN 9; skipped if the binary is not built | Phase 4 |
| 14 | `src/decifer/cpp/enumerator.cpp` | Witness mode: decided-parent check also when the parent is not in T (BFS from every node the parent may become) | Phase 3 tuning step 1 |
| 15 | `src/decifer/cpp/flashlight.cpp`, `enumerator.{h,cpp}` | Labels decided in order of depth below the mutated copy; witness mode picks the frontier arc with the smallest `targetDistance` (Gabow–Myers is correct for any arc choice) | Phase 3 tuning step 2 |
| 16 | `src/decifer/cpp/enumerator.{h,cpp}` | `targetDistance` made constraint- and label-aware: static BFS distances over genotype nodes to every (state of L, label) | Phase 3 tuning step 2 |
| 17 | `src/decifer/cpp/enumerator.{h,cpp}` | `candidates(q)`: candidate nodes propagated down the decided constraint tree; an uncovered state with a decided parent must have a candidate | Phase 3 tuning step 3 (necessary condition only, so exact) |
| 18 | `src/decifer/cpp/enumerator.{h,cpp}`, main | `route(allele)` / `shortestPath`: greedy witness construction tried before the exhaustive search (built with `add`/`report`, so the witness is a valid tree); `routed` statistic | Phase 3 tuning step 4. First version had a bug (the target check read `_owner` of BFS nodes not in T, so routing never succeeded); fixed by tracking the region of each BFS node |
| 19 | `generatestatetreesmain.cpp`, `enumerator.h`, `flashlight.{h,cpp}` | `-algorithm auto`: tree and flashlight run concurrently on each set (stop flag `setStop`/`stopped` checked by the search and by flashlight's decision loops); the first to finish is kept | Both are exact, so the output does not depend on which wins; each covers the other's weak cases |
| 20 | `src/decifer/cpp/enumerator.{h,cpp}` | `route(allele, strategy, seed)`: shortest-first, longest-first, 8 seeded random orders with random tie-breaking (`ROUTE_RANDOM_ATTEMPTS`) | Phase 3 tuning step 5 |
| 21 | `src/decifer/cpp/enumerator.{h,cpp}` | `negotiate(allele)`: negotiated-congestion routing (Dijkstra over nodes; present cost ×1.3 per round, history += overuse; ≤ 200 rounds) when all parents are decided; `findWitness` first checks the search's relaxations at the empty tree | Phase 3 tuning step 6. First version (present cost ×2, 30 rounds) never converged |
| 22 | `benchmarks/statetrees/make_cn_sets.py` | `clustered` mode: clonal state with largest copy number maxcn plus states 1–2 single-copy changes from an earlier tumor state | Uniformly random sets (e.g. {1,1;1,10;2,0;7,2}) are unlike real loci; clustered sets are a better proxy for DeCiFer inputs |
| 23 | `generatestatetreesmain.cpp` | `-symmetry`: the state trees of swap(L) are derived from those of L when both are requested (sequential: if swap(L) is already present; threaded: one of each pair is computed) | Phase 2. Verified identical: shipped (1 and 3 threads, tree and flashlight), exh3, all-combinations `-maxXY 3 -maxCN 3` |
| 24 | `generatestatetreesmain.cpp` | Fix: self-symmetric sets are not "derived" (they were computed by the legacy enumerator inside `getStateTrees`); assert before deriving | Bug found by the tests' run time. **Not compiled yet** (lemon removed) |
| 25 | `test/test_statetrees.py` | Shipped-file test over algorithm ∈ {tree, flashlight, auto} × threads ∈ {1,3} × symmetry; combinations-mode test also for auto | 19 passed (with the pre-fix binary; the symmetry+threads cases took ~110 s each because of the bug in row 24) |
| 26 | `src/decifer/cpp/solver.cpp` | `std::isnan`/`std::isinf` and `#include <cmath>` (was unqualified `isnan`/`isinf`) | Build failure with Boost 1.85 from the new `decifer-cpp` env: Boost 1.75 had pulled the global names in indirectly. Same behavior |
| 27 | `src/decifer/cpp/enumerator.{h,cpp}`, main | Mandatory cell (dominator) relaxation in witness mode (`Requirement`, `requirementPath`, `mandatoryCells`, called at the end of `reachableWitness`, so also in `findWitness`'s root check); statistic `pruned_mandatory` | User request. Requirements: each uncovered state of L with a decided parent (path from the parent's region, or from the parent's candidates) and each non-L leaf (path to an allowed child of its owner). A cell is mandatory for a requirement if blocking it disconnects the requirement's node-level BFS. Paths of requirements with different owners, or from different leaves, are cell-disjoint in every valid extension (a non-L cell has one node, with one owner; leaves have disjoint subtrees), so each requirement must avoid the mandatory cells of those; iterated to a fixed point. Sound by induction (only cells that the other requirement's real path uses are blocked), so exact |
| 28 | `src/decifer/cpp/frontier.{h,cpp}` (new), `CMakeLists.txt`, main | `FrontierEnumerator`, `-algorithm frontier`: transfer-matrix DP over the lattice (see Method: frontier DP); statistics `states`, `partials`, `max_layer`, `max_partials`, `accepted` | User request: the general method for topological obstructions (planar crossings) that no relaxation of the oracle sees |
| 29 | `benchmarks/statetrees/run_sets.py` | `check BIN REFDIR`: runs BIN once on all sets a reference run finished, splits the output into per-set blocks and compares byte for byte | The per-set process start-up dominated the gate runs |
| 30 | `src/decifer/cpp/frontier.{h,cpp}` | Static label pruning (both modes): a node with label > 0 must be reachable from an allowed mutated node (m, 1); a node not in L must reach an allowed node of L. Later refined per owner: `_toChild[o]` (reverse BFS through non-L cells from the nodes of L that may be state tree children of owner o, with allowed labels), checked for every childless non-L frontier cell | Necessary conditions only, so exact (enumeration gates re-run: identical). No measurable effect in enumeration mode; it targets the constrained mode |
| 31 | `src/decifer/cpp/frontier.{h,cpp}` | Constrained mode: `enumerate(result, constraints)` and `findWitness` (existence: one partial state tree per state, stop at the first). Decisions follow the constraints (mutation only in the decided state, decided labels, every resolved state tree parent checked; pending nodes of a component must share an allowed parent) | Exact oracle for flashlight |
| 32 | `src/decifer/cpp/enumerator.{h,cpp}` | `setGrowBudget`/`undecided`: the exhaustive witness search stops after a budget of search nodes and reports *undecided* | So the DP can take over |
| 33 | `src/decifer/cpp/flashlight.{h,cpp}`, main | Oracle → INFEASIBLE / FEASIBLE / UNDECIDED (`GROW_BUDGET` = 10,000 search nodes). UNDECIDED in the parent phase: `solve(C)` = constrained DP enumeration of every state tree consistent with the decisions so far, inserted directly, subtree skipped. UNDECIDED in the label phase: one DP for the whole label phase (parents decided, labels free). Statistics `frontier_calls`, `frontier_states`, `frontier_max_layer` | {1,1;3,5;3,6;5,5} (CN 6, never finished before): DP as existence oracle 104 s → per-decision DP enumeration 117 s → label-phase escalation **42 s**, all 471 trees, identical output across versions. Exact: the subtree's output is by definition the set of state trees consistent with its decisions |
| 34 | `src/decifer/cpp/frontier.{h,cpp}` | `initOrder()`: the processing order computed in a function. Tried and **reverted**: sweeping from the mutation's row up then down (frontier doubles: {0,0;0,4;1,1;4,4} 4.5 s → timeout) and sweeping top-down when the mutation is in the upper half (stuck set 42 s → 203 s) | Bottom-up from the root's side is best in every test |
| 35 | `src/decifer/cpp/frontier.{h,cpp}` | No allocations per transition: `Work` objects reused by recursion depth (`_pool`), keys stored once (in the layer's hash map), one reused key buffer, partial state trees of a layer in one flat vector, sorted and deduplicated per step | `perf`: malloc/free were ~40% of the time. {1,1;5,5} 5.9 → 3.6 s; CN-6 stuck set 42 → 25 s; output identical |
| 36 | `src/decifer/cpp/frontier.{h,cpp}` | `HAS_CHILD` dropped from the key for nodes of L and the mutation cell; open-end deficit bound: (# childless non-L frontier nodes) − (# components without the root) ≤ (# unprocessed states of L) + (# pending nodes) | Exact: open ends of a component have disjoint subtrees, so each needs its own first descendant in L, and two can share one only if one's component attaches below the other (once per component). {1,1;5,5} 2.57 → 1.64 M states (1.9 s); all gates identical |
| 37 | `frontier.{h,cpp}`, `flashlight.{h,cpp}`, main | State limit per DP step (`setStateLimit`/`exceeded`, option `-frontierLimit`, default 2,000,000 ≈ 1 GB at the measured ~0.5 kB per state of the largest step). Flashlight: a DP over the limit → the decision is made by the unlimited exhaustive search (`decided`, `_unbounded` for the label phase); statistic `frontier_exceeded` | A CN-8 constrained DP ({1,1;6,7;6,8;7,8}, width 7) keeps ~640 k states per step for minutes; wider lattices would exhaust memory |
| 38 | main | `auto` races tree, flashlight and frontier (DP with the state limit; it gives up rather than win when over it); 3 threads per set | frontier is the fastest engine on narrow lattices: random thin sets 0.5 s (tree 21.9 s), random square 2.1 s (tree 58.9 s) |
| 39 | `test/test_statetrees.py` | `frontier` in the shipped/legacy/combinations/beyond-shipped tests; new `test_flashlight_dynamic_program_fallback` (two exh4 sets where flashlight calls the DP, with `-frontierLimit` 0, 1 and 2,000,000; 1 forces the over-limit fallback) | 27 passed |
| 40 | main | Default algorithm `tree` → `auto` | Byte-identical everywhere; at worst ~3× slower than the fastest engine (on one core), while a single engine can be orders of magnitude slower |
| 41 | `scripts/README.md` | Algorithms, `-threads`, `-frontierLimit`, `-stats`; measured run times replace "if `MAX_CN` is greater than 6, `generatestatetrees` may take a very long time" | Phase 4 (user-facing docs) |
| 42 | `benchmarks/statetrees/sets/` (new) | `clustered.txt`, `rand_square.txt`, `rand_thin.txt` verbatim, with a README | Reproducibility: the commands/seeds that made them were not recorded and the current generator does not reproduce them with seeds 0–300 (checked); the other lists (shipped, exhaustive, ladder) are reproduced exactly by `make_cn_sets.py` (checked) |
| 43 | `src/decifer/cpp/enumerator.h` | Class comment: `findWitness` (not the nonexistent `setConstraints`), the witness-mode relaxations and the budget | Accuracy |
| — | (incident) | I rebuilt `$TMPDIR/st-dev/generatestatetrees` (allocation removal, row 35) while the clustered benchmark `dp-flashlight-clustered` was ending, so its last sets may have used the new build (same behavior, different speed). Benchmarks since use frozen copies (`gst-dp2`, `gst-dp3`, `gst-dp4`) | Transparency |
| — | (test) | `$TMPDIR/st-b1/gst`: copy with `GROW_BUDGET = 1`, so every oracle call that routing does not decide goes to the constrained DP: shipped 469/469 and exh3 576/576 identical, with 4,011 and 4,946 DP calls | Validates constrained mode and its integration |
| — | (incident) | A test loop ran `-algorithm tree` without a timeout on {1,1;2,2;7,6}; I stopped the background command after > 10 min. The sandbox hides processes of other commands, so I could not confirm that the orphaned `generatestatetrees` child was killed | Transparency |
| — | (incident) | I rebuilt `$TMPDIR/st-new/generatestatetrees` while a background run (flashlight ladder + exh4) was using it, so later sets in that run used the newer binary. Only main's dispatch had changed. All development builds since go to `$TMPDIR/st-dev`, and the gate will be rerun on the final binary | Transparency |

## Sandbox-limited tasks and workarounds

| Task | Blocked by | Workaround |
|---|---|---|
| Rebuilding `generatestatetrees` after 16:44 on 2026-10-21 | The conda package cache `/home/myersm2/miniforge3/pkgs` was cleaned at 16:44:52 (directory mtime). `lemon-1.3.1-h9457a4f_11` (the only lemon headers on the machine), `libboost-devel-1.85.0` and `cmake-3.31.8` are gone. **Not done by me.** No network to re-download | **Resolved 2026-10-22 by the user**: env `/home/myersm2/miniforge3/envs/decifer-cpp` (lemon 1.3.1, Boost 1.85, cmake, make), outside the package cache so `conda clean` cannot remove it. With it, `solver.cpp` needed `std::isnan`/`std::isinf` (row 26) |
| Recovering the commands that made the clustered/random set lists from the previous session's transcript (`/home/myersm2/.claude/projects/-data1-shahs3-users-myersm2-repos-decifer/25c2bbb0-46f7-4b56-b43f-136576a939de.jsonl`) | Not readable from the sandbox (reported as nonexistent) | Saved the exact lists in `benchmarks/statetrees/sets/` instead (row 42) |
| (observation) `git status` inside my sandbox lists untracked `.bashrc`, `.gitconfig`, `.mcp.json`, `.vscode`, `.claude/` (and similar) in the repository root | They are the sandbox's protective mounts: character devices 1,3 (`/dev/null`) dated 2026-10-15, plus `.claude/` holding only such stubs and the harness's `.cc-writes`. Not created or edited by me | Left untouched; do not commit them |
| `mamba create -n decifer-cpp --offline cmake make gxx libboost-devel lemon` | No network (proxy refuses `conda.anaconda.org`); the solver picked uncached versions | Pinned to cached versions; still failed because `liblzma-devel` is not cached. Removed the empty env stub the failed attempt left behind (`envs/decifer-cpp`, contained only `conda-meta/history`). Used system cmake/g++ with cached lemon headers and the `hatchet` env's Boost (read-only) |
