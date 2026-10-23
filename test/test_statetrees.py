"""
Tests of the C++ state tree enumeration (`generatestatetrees`).

The binary is looked up in $GENERATESTATETREES, then build/generatestatetrees (the build directory
of the README instructions), then on PATH. The tests are skipped if it is not found.

All algorithms (`legacy`, the original enumeration; `tree`; `flashlight`; `frontier`; `auto`, which
runs tree, flashlight and frontier concurrently) and options (`-threads`, `-symmetry`,
`-frontierLimit`) must give byte-identical output, and must reproduce the shipped state_trees.txt.
"""

import os
import shutil
import subprocess

import pytest


this_dir = os.path.dirname(__file__)
REPO = os.path.dirname(this_dir)
SHIPPED = os.path.join(REPO, "src", "decifer", "state_trees.txt")


def find_binary():
    candidates = [os.environ.get("GENERATESTATETREES"),
                  os.path.join(REPO, "build", "generatestatetrees"),
                  shutil.which("generatestatetrees")]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


BINARY = find_binary()
pytestmark = pytest.mark.skipif(BINARY is None, reason="generatestatetrees is not built")


def shipped_sets():
    """Copy number state sets of the shipped state tree file, in file order."""
    with open(SHIPPED) as f:
        lines = f.read().split("\n")
    sets, i = [], 1
    for _ in range(int(lines[0].split()[0])):
        sets.append(lines[i].split())
        i += 2 + int(lines[i + 1].split()[0])
    return sets


def run(tmp_path, sets, *args):
    """Run generatestatetrees on the given sets (lists of "x,y" strings) and return its stdout."""
    inp = tmp_path / "cn_states.txt"
    inp.write_text("".join(";".join(s) + "\n" for s in sets))
    out = subprocess.run([BINARY] + list(args) + [str(inp)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         check=True)
    return out.stdout.decode()


@pytest.mark.parametrize("algorithm", ["tree", "flashlight", "frontier", "auto"])
@pytest.mark.parametrize("threads", ["1", "3"])
@pytest.mark.parametrize("symmetry", [[], ["-symmetry"]])
def test_reproduces_shipped_file(tmp_path, algorithm, threads, symmetry):
    with open(SHIPPED) as f:
        expected = f.read()
    assert run(tmp_path, shipped_sets(), "-algorithm", algorithm, "-threads", threads, *symmetry) == expected


def test_algorithms_agree_with_legacy(tmp_path):
    # sets that the original enumeration finishes quickly (maximum copy number at most 2)
    sets = [s for s in shipped_sets() if max(int(c) for xy in s for c in xy.split(",")) <= 2]
    assert len(sets) > 50
    legacy = run(tmp_path, sets, "-algorithm", "legacy")
    assert run(tmp_path, sets, "-algorithm", "tree") == legacy
    assert run(tmp_path, sets, "-algorithm", "flashlight") == legacy
    assert run(tmp_path, sets, "-algorithm", "frontier") == legacy


@pytest.mark.parametrize("algorithm", ["tree", "flashlight", "frontier", "auto"])
def test_all_combinations_mode_agrees_with_legacy(algorithm):
    # every pair of copy number states with at most 2 copies, including sets without (1,1)
    def combos(*args):
        return subprocess.run([BINARY, "-maxXY", "2", "-maxCN", "2"] + list(args), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, check=True).stdout.decode()
    assert combos("-algorithm", algorithm) == combos("-algorithm", "legacy")


def test_known_state_trees(tmp_path):
    out = run(tmp_path, [["1,1", "3,3"]], "-algorithm", "flashlight")
    assert out == ("1 #copynumber_states\n"
                   "1,1 3,3\n"
                   "3 #state_trees\n"
                   "1,1,0 1,1,1 1,1,0 3,3,0\n"
                   "1,1,0 1,1,1 1,1,1 3,3,3\n"
                   "1,1,0 3,3,0 3,3,0 3,3,1\n")
    assert run(tmp_path, [["0,0", "1,0", "1,1"]]).split("\n")[2] == "12 #state_trees"
    assert run(tmp_path, [["1,1", "2,2", "2,3", "3,2"]]).split("\n")[2] == "228 #state_trees"


def test_tree_and_flashlight_agree_beyond_shipped_file(tmp_path):
    # copy numbers up to 5, where the original enumeration takes minutes to hours
    sets = [["1,1", "5,5"], ["1,1", "2,2", "5,4"], ["0,0", "0,4", "1,1", "4,0"], ["1,1", "2,2", "2,4", "4,2"],
            ["1,1", "5,0", "5,2"], ["1,1", "3,1", "5,1"]]
    tree = run(tmp_path, sets, "-algorithm", "tree")
    assert run(tmp_path, sets, "-algorithm", "flashlight") == tree
    assert run(tmp_path, sets, "-algorithm", "frontier") == tree


@pytest.mark.parametrize("limit", ["0", "1", "2000000"])
def test_flashlight_dynamic_program_fallback(tmp_path, limit):
    # sets where the exhaustive search of flashlight's oracle exceeds its budget, so the frontier
    # dynamic program enumerates the state trees of the undecided prefix; with a state limit of 1
    # it gives up and the unlimited exhaustive search decides instead
    sets = [["1,1", "1,4", "2,0", "4,2"], ["1,1", "3,0", "4,0", "4,4"]]
    assert run(tmp_path, sets, "-algorithm", "flashlight", "-frontierLimit", limit) == \
        run(tmp_path, sets, "-algorithm", "tree")


def test_flashlight_high_copy_numbers(tmp_path):
    # far beyond the reach of the subtree enumerations (both alleles amplified to 9)
    out = run(tmp_path, [["1,1", "9,9"]], "-algorithm", "flashlight")
    assert out.split("\n")[2] == "3 #state_trees"
    assert "1,1,1 9,9,9" in out
