/*
 * flashlight.cpp
 *
 *  Enumeration of state trees by backtracking over the decisions that define a state tree,
 *  with a feasibility oracle (flashlight search).
 */

#include "flashlight.h"

FlashlightEnumerator::FlashlightEnumerator(const IntPairSet& L, uint64_t frontierLimit)
  : _E(L)
  , _D(L)
  , _nL(L.size())
  , _rootIndex(-1)
  , _labelNodes()
  , _unbounded(false)
  , _result(NULL)
  , _stats()
{
  for (int i = 0; i < _nL; ++i)
  {
    if (_E.states()[i] == IntPair(1, 1))
    {
      _rootIndex = i;
    }
  }
  _D.setStateLimit(frontierLimit);
}

FlashlightEnumerator::OracleResult FlashlightEnumerator::oracle(const Constraints& C, Witness& W, bool unbounded)
{
  ++_stats._oracleCalls;
  _E.setGrowBudget(unbounded || _unbounded ? 0 : GROW_BUDGET);
  if (_E.findWitness(C, W._S, W._parents, W._labels, W._mutationState))
  {
    ++_stats._feasible;
    return FEASIBLE;
  }
  return _E.undecided() ? UNDECIDED : INFEASIBLE;
}

bool FlashlightEnumerator::solve(const Constraints& C)
{
  ++_stats._frontierCalls;
  StateGraph::StateEdgeSetSet trees;
  _D.enumerate(trees, &C);
  if (_D.exceeded())
  {
    ++_stats._frontierExceeded;
    return false;
  }
  _result->insert(trees.begin(), trees.end());
  return true;
}

void FlashlightEnumerator::decided(OracleResult r, int i, Constraints& C, Witness& W)
{
  if (r == UNDECIDED && !solve(C))
  {
    // the dynamic program needs too much memory: decide by the unlimited exhaustive search
    r = oracle(C, W, true);
  }
  if (r == FEASIBLE)
  {
    decideParent(i, C, W);
  }
}

bool FlashlightEnumerator::cycle(const Constraints& C, int i) const
{
  // any other cycle would have been rejected when its last parent was decided
  int j = C._parent[i];
  for (int steps = 0; j >= 0 && steps <= _nL + 1; ++steps)
  {
    if (j == i)
    {
      return true;
    }
    j = C._parent[j];
  }
  return false;
}

void FlashlightEnumerator::enumerate(StateGraph::StateEdgeSetSet& result)
{
  _result = &result;

  for (int m = 0; m < _nL && !_E.stopped(); ++m)
  {
    const IntPair& xy = _E.states()[m];
    if (xy.first == 0 && xy.second == 0) continue;

    // the mutated copy (index nL) is a child of m with label 1, the root (1,1) has no parent,
    // and neither the root nor m (before the mutation) carry the mutation
    Constraints C;
    C._mutationState = m;
    C._parent = IntVector(_nL + 1, -1);
    C._label = IntVector(_nL + 1, -1);
    C._parent[_nL] = m;
    C._label[_nL] = 1;
    C._label[m] = 0;
    if (_rootIndex != -1)
    {
      C._parent[_rootIndex] = -2;
      C._label[_rootIndex] = 0;
    }

    Witness W;
    decided(oracle(C, W), 0, C, W);
  }

  _result = NULL;
}

void FlashlightEnumerator::decideParent(int i, Constraints& C, const Witness& W)
{
  if (_E.stopped())
  {
    return;
  }
  
  if (i == _nL)
  {
    // labels: decided for descendants of the mutated copy, 0 for all other nodes
    const IntVector labels = C._label;
    // decided parents first: order the descendants by their depth below the mutated copy
    std::vector<IntPair> depthNode;
    for (int j = 0; j < _nL; ++j)
    {
      int k = C._parent[j];
      int depth = 1;
      while (k >= 0 && k != _nL)
      {
        k = C._parent[k];
        ++depth;
      }
      if (k == _nL)
      {
        depthNode.push_back(IntPair(depth, j));
      }
      else
      {
        C._label[j] = 0;
      }
    }
    std::sort(depthNode.begin(), depthNode.end());
    _labelNodes.clear();
    for (const IntPair& dj : depthNode)
    {
      _labelNodes.push_back(dj.second);
    }
    if (decideLabel(0, C, W) && !solve(C))
    {
      // all labels at once by the dynamic program needs too much memory: decide them by the
      // unlimited exhaustive search
      _unbounded = true;
      decideLabel(0, C, W);
      _unbounded = false;
    }
    C._label = labels;
    return;
  }

  if (C._parent[i] != -1)
  {
    // the root
    decideParent(i + 1, C, W);
    return;
  }

  // -2: no parent, only possible for nodes without an ancestor in L, i.e. if (1,1) is not in L
  for (int p = -2; p <= _nL; ++p)
  {
    if (p == -1 || p == i || (p == -2 && _rootIndex != -1)) continue;

    C._parent[i] = p;
    if (cycle(C, i))
    {
      ++_stats._cycles;
      continue;
    }

    if (W._parents[i] == p)
    {
      ++_stats._reused;
      decideParent(i + 1, C, W);
    }
    else
    {
      Witness W2;
      decided(oracle(C, W2), i + 1, C, W2);
    }
  }
  C._parent[i] = -1;
}

bool FlashlightEnumerator::decideLabel(int k, Constraints& C, const Witness& W)
{
  if (_E.stopped())
  {
    return false;
  }
  
  if (k == static_cast<int>(_labelNodes.size()))
  {
    // every decision is made, so the witness's state tree is the one they define
    _result->insert(W._S);
    return false;
  }

  const int j = _labelNodes[k];
  const IntPair& xy = _E.states()[j];
  const int maxLabel = std::max(xy.first, xy.second);
  for (int z = 0; z <= maxLabel; ++z)
  {
    C._label[j] = z;
    bool undecided = false;
    if (W._labels[j] == z)
    {
      ++_stats._reused;
      undecided = decideLabel(k + 1, C, W);
    }
    else
    {
      Witness W2;
      const OracleResult r = oracle(C, W2);
      undecided = r == UNDECIDED || (r == FEASIBLE && decideLabel(k + 1, C, W2));
    }
    if (undecided)
    {
      C._label[j] = -1;
      return true;
    }
  }
  C._label[j] = -1;
  return false;
}
