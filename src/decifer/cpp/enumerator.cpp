/*
 * enumerator.cpp
 *
 *  Exact state tree enumeration with pruning of partial trees that cannot become valid.
 */

#include "enumerator.h"
#include <limits>
#include <algorithm>
#include <queue>
#include <functional>

StateTreeEnumerator::StateTreeEnumerator(const IntPairSet& L)
  : _L(L)
  , _states(L.begin(), L.end())
  , _nL(L.size())
  , _maxX(1)
  , _maxY(1)
  , _nCells(0)
  , _root(-1)
  , _nonLLeaves(0)
  , _coveredL(0)
  , _mutationArc(-1)
  , _stamp(0)
  , _stampNode(0)
  , _witnessMode(false)
  , _found(false)
  , _witnessMutationState(-1)
  , _result(NULL)
  , _stop(NULL)
  , _growBudget(0)
  , _growCalls(0)
  , _undecided(false)
{
  // same state graph dimensions as StateGraph::getStateTrees
  for (const IntPair& xy : L)
  {
    _maxX = std::max(_maxX, xy.first);
    _maxY = std::max(_maxY, xy.second);
  }
  _nCells = (_maxX + 1) * (_maxY + 1);

  _lIndex = IntVector(_nCells, -1);
  for (int i = 0; i < _nL; ++i)
  {
    _lIndex[cell(_states[i].first, _states[i].second)] = i;
  }

  // copy number moves: an amplification of a copy requires 0 < x < max_x, a deletion x > 0.
  // For every node (x, y, xbar, ybar) one of the two variants (mutated or nonmutated copy) of
  // each move exists, so reachability between copy number states does not depend on xbar, ybar.
  _cellNbrs = IntMatrix(_nCells);
  for (int x = 0; x <= _maxX; ++x)
  {
    for (int y = 0; y <= _maxY; ++y)
    {
      IntVector& nbrs = _cellNbrs[cell(x, y)];
      if (0 < x && x < _maxX) nbrs.push_back(cell(x + 1, y));
      if (0 < y && y < _maxY) nbrs.push_back(cell(x, y + 1));
      if (x > 0) nbrs.push_back(cell(x - 1, y));
      if (y > 0) nbrs.push_back(cell(x, y - 1));
    }
  }

  // copy number states from which some state of L is reachable
  std::vector<char> reachesL(_nCells, 0);
  bool changed = true;
  for (int c = 0; c < _nCells; ++c)
  {
    reachesL[c] = _lIndex[c] >= 0;
  }
  while (changed)
  {
    changed = false;
    for (int c = 0; c < _nCells; ++c)
    {
      if (reachesL[c]) continue;
      for (int d : _cellNbrs[c])
      {
        if (reachesL[d])
        {
          reachesL[c] = 1;
          changed = true;
          break;
        }
      }
    }
  }

  // nodes, as in StateGraph::init, except those that cannot reach L
  _toNode = IntVector(_nCells * (_maxX + 1) * (_maxY + 1), -1);
  for (int x = 0; x <= _maxX; ++x)
  {
    for (int y = 0; y <= _maxY; ++y)
    {
      if (!reachesL[cell(x, y)]) continue;
      for (int xbar = 0; xbar <= x; ++xbar)
      {
        _toNode[cell(x, y) * (_maxX + 1) * (_maxY + 1) + xbar * (_maxY + 1)] = _x.size();
        _x.push_back(x);
        _y.push_back(y);
        _xbar.push_back(xbar);
        _ybar.push_back(0);
        _cell.push_back(cell(x, y));
      }
      for (int ybar = 1; ybar <= y; ++ybar)
      {
        _toNode[cell(x, y) * (_maxX + 1) * (_maxY + 1) + ybar] = _x.size();
        _x.push_back(x);
        _y.push_back(y);
        _xbar.push_back(0);
        _ybar.push_back(ybar);
        _cell.push_back(cell(x, y));
      }
    }
  }
  _root = _toNode[cell(1, 1) * (_maxX + 1) * (_maxY + 1)];
  assert(_root != -1);

  // arcs, as in StateGraph::init, except mutation arcs from copy number states not in L
  const int nNodes = _x.size();
  _outArcs = IntMatrix(nNodes);
  for (int v = 0; v < nNodes; ++v)
  {
    const int x = _x[v];
    const int y = _y[v];
    const int xbar = _xbar[v];
    const int ybar = _ybar[v];

    // mutation edges
    if (xbar == 0 && ybar == 0 && x > 0 && _lIndex[_cell[v]] >= 0)
    {
      addArc(v, MUTATION, x, y, xbar + 1, ybar);
    }
    if (xbar == 0 && ybar == 0 && y > 0 && _lIndex[_cell[v]] >= 0)
    {
      addArc(v, MUTATION, x, y, xbar, ybar + 1);
    }
    // amplification edges
    if (0 < x && x < _maxX && xbar < x)
    {
      addArc(v, AMPLIFICATION, x+1, y, xbar, ybar);
    }
    if (0 < x && x < _maxX && xbar > 0)
    {
      addArc(v, AMPLIFICATION, x+1, y, xbar+1, ybar);
    }
    if (0 < y && y < _maxY && ybar < y)
    {
      addArc(v, AMPLIFICATION, x, y+1, xbar, ybar);
    }
    if (0 < y && y < _maxY && ybar > 0)
    {
      addArc(v, AMPLIFICATION, x, y+1, xbar, ybar+1);
    }
    // deletion edges
    if (x > xbar)
    {
      addArc(v, DELETION, x-1, y, xbar, ybar);
    }
    if (xbar > 0)
    {
      addArc(v, DELETION, x-1, y, xbar-1, ybar);
    }
    if (y > ybar)
    {
      addArc(v, DELETION, x, y-1, xbar, ybar);
    }
    if (ybar > 0)
    {
      addArc(v, DELETION, x, y-1, xbar, ybar-1);
    }
  }

  // number of arcs from each node to a node of each state of L with each label (search
  // guidance in witness mode), by breadth-first search over the reversed arcs
  IntMatrix inArcs(nNodes);
  for (int a = 0; a < static_cast<int>(_arcs.size()); ++a)
  {
    inArcs[_arcs[a]._target].push_back(a);
  }
  _targetOffset = IntVector(_nL + 1, 0);
  for (int i = 0; i < _nL; ++i)
  {
    _targetOffset[i + 1] = _targetOffset[i] + std::max(_states[i].first, _states[i].second) + 1;
  }
  const int unreachable = std::numeric_limits<int>::max();
  _distToTarget = IntMatrix(_targetOffset[_nL], IntVector(nNodes, unreachable));
  for (int i = 0; i < _nL; ++i)
  {
    const int c = cell(_states[i].first, _states[i].second);
    for (int z = 0; z <= std::max(_states[i].first, _states[i].second); ++z)
    {
      IntVector& dist = _distToTarget[_targetOffset[i] + z];
      IntVector queue;
      for (int v = 0; v < nNodes; ++v)
      {
        if (_cell[v] == c && std::max(_xbar[v], _ybar[v]) == z)
        {
          dist[v] = 0;
          queue.push_back(v);
        }
      }
      for (size_t h = 0; h < queue.size(); ++h)
      {
        for (int a : inArcs[queue[h]])
        {
          const int u = _arcs[a]._source;
          if (dist[u] == unreachable)
          {
            dist[u] = dist[queue[h]] + 1;
            queue.push_back(u);
          }
        }
      }
    }
  }
  
  _cellNodes = IntMatrix(_nCells);
  for (int v = 0; v < nNodes; ++v)
  {
    _cellNodes[_cell[v]].push_back(v);
  }
  
  _F = IntMatrix(nNodes + 2);
  _FF = IntMatrix(nNodes + 2);
  _seen = IntVector(_nCells, 0);
  _seenNode = IntVector(nNodes, 0);
  _predArc = IntVector(nNodes, -1);
  _bfsOwnerIndex = IntVector(nNodes, -2);
  _blocked = std::vector<char>(_nCells, 0);
  _cand = IntMatrix(_nL + 1);
  _candState = IntVector(_nL + 1, 0);
}

void StateTreeEnumerator::addArc(int v, ArcType type, int x, int y, int xbar, int ybar)
{
  int w = _toNode[cell(x, y) * (_maxX + 1) * (_maxY + 1) + xbar * (_maxY + 1) + ybar];
  if (w == -1) return;

  StateArc a;
  a._source = v;
  a._target = w;
  a._type = type;
  _outArcs[v].push_back(_arcs.size());
  _arcs.push_back(a);
}

void StateTreeEnumerator::init()
{
  const int nNodes = _x.size();
  _inT = std::vector<char>(nNodes, 0);
  _disabled = std::vector<char>(_arcs.size(), 0);
  _parent = IntVector(nNodes, -1);
  _owner = IntVector(nNodes, -1);
  _childCount = IntVector(nNodes, 0);
  _cellUsed = IntVector(_nCells, 0);
  _tNodes.clear();
  _mutationArc = -1;
  _found = false;

  _inT[_root] = 1;
  _cellUsed[_cell[_root]] = 1;
  _owner[_root] = isL(_root) ? _root : -1;
  _coveredL = isL(_root) ? 1 : 0;
  _nonLLeaves = isL(_root) ? 0 : 1;
  _tNodes.push_back(_root);
  _nodeOfIndex = IntVector(_nL + 1, -1);
  if (isL(_root))
  {
    _nodeOfIndex[_lIndex[_cell[_root]]] = _root;
  }

  _F[0] = _outArcs[_root];
}

void StateTreeEnumerator::add(int uv)
{
  const int u = _arcs[uv]._source;
  const int v = _arcs[uv]._target;

  _inT[v] = 1;
  _parent[v] = u;
  if (_childCount[u]++ == 0 && !isL(u))
  {
    --_nonLLeaves;
  }
  if (!isL(v))
  {
    ++_nonLLeaves;
  }
  if (_cellUsed[_cell[v]]++ == 0 && isL(v))
  {
    ++_coveredL;
  }
  _owner[v] = isL(v) ? v : _owner[u];
  if (_arcs[uv]._type == MUTATION)
  {
    _mutationArc = uv;
  }
  if (isL(v))
  {
    _nodeOfIndex[nodeIndex(v)] = v;
  }
  _tNodes.push_back(v);
}

void StateTreeEnumerator::remove(int uv)
{
  const int u = _arcs[uv]._source;
  const int v = _arcs[uv]._target;

  assert(_tNodes.back() == v);
  _tNodes.pop_back();
  if (isL(v))
  {
    _nodeOfIndex[nodeIndex(v)] = -1;
  }
  if (uv == _mutationArc)
  {
    _mutationArc = -1;
  }
  _owner[v] = -1;
  if (--_cellUsed[_cell[v]] == 0 && isL(v))
  {
    --_coveredL;
  }
  if (!isL(v))
  {
    --_nonLLeaves;
  }
  if (--_childCount[u] == 0 && !isL(u))
  {
    ++_nonLLeaves;
  }
  _parent[v] = -1;
  _inT[v] = 0;
}

bool StateTreeEnumerator::mutationAvailable(const IntVector& F) const
{
  // an uncovered state of L where the mutation may occur
  for (int i = 0; i < _nL; ++i)
  {
    if (_witnessMode && _constraints._mutationState >= 0 && _constraints._mutationState != i) continue;
    const IntPair& xy = _states[i];
    if ((xy.first > 0 || xy.second > 0) && _cellUsed[cell(xy.first, xy.second)] == 0)
    {
      return true;
    }
  }
  // or a mutation arc in the frontier (its source is in L)
  for (int a : F)
  {
    if (_arcs[a]._type != MUTATION) continue;
    if (_witnessMode && _constraints._mutationState >= 0
        && _lIndex[_cell[_arcs[a]._source]] != _constraints._mutationState) continue;
    return true;
  }
  return false;
}

bool StateTreeEnumerator::reachable()
{
  if (_witnessMode)
  {
    return reachableWitness();
  }
  
  // every uncovered state of L must be reachable from T through unused copy number states
  ++_stamp;
  _queue.clear();
  for (int v : _tNodes)
  {
    const int c = _cell[v];
    if (_seen[c] != _stamp)
    {
      _seen[c] = _stamp;
      _queue.push_back(c);
    }
  }
  for (size_t h = 0; h < _queue.size(); ++h)
  {
    for (int d : _cellNbrs[_queue[h]])
    {
      if (_cellUsed[d] == 0 && _seen[d] != _stamp)
      {
        _seen[d] = _stamp;
        _queue.push_back(d);
      }
    }
  }
  for (int i = 0; i < _nL; ++i)
  {
    const int c = cell(_states[i].first, _states[i].second);
    if (_cellUsed[c] == 0 && _seen[c] != _stamp)
    {
      return false;
    }
  }

  // every leaf not in L must reach, through unused copy number states not in L, an uncovered
  // state of L that may become a child (in the state tree) of the leaf's nearest ancestor in L
  for (int v : _tNodes)
  {
    if (_childCount[v] > 0 || isL(v)) continue;

    int ownerIndex = -2;
    if (_witnessMode && _owner[v] != -1)
    {
      ownerIndex = nodeIndex(_owner[v]);
    }

    ++_stamp;
    _queue.clear();
    _seen[_cell[v]] = _stamp;
    _queue.push_back(_cell[v]);
    bool found = false;
    for (size_t h = 0; h < _queue.size() && !found; ++h)
    {
      for (int d : _cellNbrs[_queue[h]])
      {
        if (_cellUsed[d] != 0 || _seen[d] == _stamp) continue;
        _seen[d] = _stamp;
        const int i = _lIndex[d];
        if (i >= 0)
        {
          if (!_witnessMode || _constraints._parent[i] == -1 || _constraints._parent[i] == ownerIndex)
          {
            found = true;
            break;
          }
        }
        else
        {
          _queue.push_back(d);
        }
      }
    }
    if (!found)
    {
      return false;
    }
  }

  return true;
}

bool StateTreeEnumerator::reachesChild(int target, int ownerIndex)
{
  for (size_t h = 0; h < _queueNode.size(); ++h)
  {
    const int u = _queueNode[h];
    for (int uw : _outArcs[u])
    {
      if (_arcs[uw]._type == MUTATION) continue;
      const int w = _arcs[uw]._target;
      const int cw = _cell[w];
      if (_cellUsed[cw] != 0 || _seenNode[w] == _stampNode) continue;
      _seenNode[w] = _stampNode;
      const int i = _lIndex[cw];
      if (i < 0)
      {
        _queueNode.push_back(w);
      }
      else if ((target == -1 || target == i)
               && (_constraints._parent[i] == -1 || _constraints._parent[i] == ownerIndex)
               && (_constraints._label[i] == -1 || _constraints._label[i] == std::max(_xbar[w], _ybar[w])))
      {
        return true;
      }
    }
  }
  return false;
}

bool StateTreeEnumerator::reachableWitness()
{
  // every uncovered state of L must be reachable from T through unused copy number states
  ++_stamp;
  _queue.clear();
  for (int v : _tNodes)
  {
    const int c = _cell[v];
    if (_seen[c] != _stamp)
    {
      _seen[c] = _stamp;
      _queue.push_back(c);
    }
  }
  for (size_t h = 0; h < _queue.size(); ++h)
  {
    for (int d : _cellNbrs[_queue[h]])
    {
      if (_cellUsed[d] == 0 && _seen[d] != _stamp)
      {
        _seen[d] = _stamp;
        _queue.push_back(d);
      }
    }
  }
  for (int i = 0; i < _nL; ++i)
  {
    const int c = cell(_states[i].first, _states[i].second);
    if (_cellUsed[c] == 0 && _seen[c] != _stamp)
    {
      return false;
    }
  }
  
  // every leaf not in L must reach, through nodes in unused copy number states not in L, a node
  // of an uncovered state of L whose parent and label constraints allow it to be a state tree
  // child of the leaf's nearest ancestor in L
  for (int v : _tNodes)
  {
    if (_childCount[v] > 0 || isL(v)) continue;
    
    ++_stampNode;
    _queueNode.clear();
    _seenNode[v] = _stampNode;
    _queueNode.push_back(v);
    if (!reachesChild(-1, _owner[v] == -1 ? -2 : nodeIndex(_owner[v])))
    {
      return false;
    }
  }
  
  // every uncovered state of L with a required state tree parent must have a candidate node
  _candState.assign(_nL + 1, 0);
  for (int i = 0; i < _nL; ++i)
  {
    if (_constraints._parent[i] == -1 || _cellUsed[cell(_states[i].first, _states[i].second)] != 0) continue;
    if (candidates(i).empty())
    {
      return false;
    }
  }

  if (!mandatoryCells())
  {
    ++_stats._prunedMandatory;
    return false;
  }

  return true;
}

bool StateTreeEnumerator::requirementPath(const Requirement& r, IntVector* cells)
{
  ++_stampNode;
  _queueNode.clear();
  for (int v : r._sources)
  {
    if (_seenNode[v] == _stampNode) continue;
    _seenNode[v] = _stampNode;
    _predArc[v] = -1;
    _queueNode.push_back(v);
  }

  int found = -1;
  for (size_t h = 0; h < _queueNode.size() && found == -1; ++h)
  {
    for (int uw : _outArcs[_queueNode[h]])
    {
      if (_arcs[uw]._type == MUTATION) continue;
      const int w = _arcs[uw]._target;
      const int cw = _cell[w];
      if (_cellUsed[cw] != 0 || _blocked[cw] || _seenNode[w] == _stampNode) continue;
      _seenNode[w] = _stampNode;
      _predArc[w] = uw;
      const int i = _lIndex[cw];
      if (i < 0)
      {
        _queueNode.push_back(w);
      }
      else if ((r._target == -1
                ? _constraints._parent[i] == -1 || _constraints._parent[i] == r._owner
                : i == r._target)
               && labelAllowed(i, w))
      {
        found = w;
        break;
      }
    }
  }
  if (found == -1)
  {
    return false;
  }

  if (cells != NULL)
  {
    cells->clear();
    for (int w = _arcs[_predArc[found]]._source; _predArc[w] != -1; w = _arcs[_predArc[w]]._source)
    {
      cells->push_back(_cell[w]);
    }
  }
  return true;
}

bool StateTreeEnumerator::mandatoryCells()
{
  // requirements: every uncovered state of L with a decided parent, from the parent's region (or
  // from the parent's candidates if the parent is not in T), and every leaf not in L
  _requirements.clear();
  for (int i = 0; i < _nL; ++i)
  {
    const int q = _constraints._parent[i];
    if (q == -1 || _cellUsed[cell(_states[i].first, _states[i].second)] != 0) continue;
    Requirement r;
    r._owner = q;
    r._leaf = -1;
    r._target = i;
    if (q >= 0 && _nodeOfIndex[q] == -1)
    {
      r._sources = candidates(q);
    }
    else
    {
      const int qNode = q == -2 ? -1 : _nodeOfIndex[q];
      for (int v : _tNodes)
      {
        if (_owner[v] == qNode)
        {
          r._sources.push_back(v);
        }
      }
    }
    _requirements.push_back(r);
  }
  for (int v : _tNodes)
  {
    if (_childCount[v] > 0 || isL(v)) continue;
    Requirement r;
    r._owner = _owner[v] == -1 ? -2 : nodeIndex(_owner[v]);
    r._leaf = v;
    r._target = -1;
    r._sources.push_back(v);
    _requirements.push_back(r);
  }

  // a single requirement is covered by the other relaxations
  const int n = _requirements.size();
  if (n < 2)
  {
    return true;
  }

  // mandatory cells of each requirement, avoiding the mandatory cells of the requirements whose
  // paths must be disjoint from it; a requirement is recomputed when those grow
  std::vector<char> dirty(n, 1);
  IntVector path;
  bool changed = true;
  while (changed)
  {
    changed = false;
    for (int k = 0; k < n; ++k)
    {
      if (!dirty[k]) continue;
      dirty[k] = 0;
      Requirement& r = _requirements[k];
      for (int l = 0; l < n; ++l)
      {
        if (l == k || !disjoint(r, _requirements[l])) continue;
        for (int c : _requirements[l]._mandatory)
        {
          _blocked[c] = 1;
        }
      }

      const bool ok = requirementPath(r, &path);
      bool grown = false;
      if (ok)
      {
        for (int c : path)
        {
          if (std::find(r._mandatory.begin(), r._mandatory.end(), c) != r._mandatory.end()) continue;
          _blocked[c] = 1;
          if (!requirementPath(r, NULL))
          {
            r._mandatory.push_back(c);
            grown = true;
          }
          _blocked[c] = 0;
        }
      }

      for (int l = 0; l < n; ++l)
      {
        for (int c : _requirements[l]._mandatory)
        {
          _blocked[c] = 0;
        }
      }
      if (!ok)
      {
        return false;
      }
      if (grown)
      {
        changed = true;
        for (int l = 0; l < n; ++l)
        {
          if (l != k && disjoint(r, _requirements[l]))
          {
            dirty[l] = 1;
          }
        }
      }
    }
  }

  return true;
}

bool StateTreeEnumerator::labelAllowed(int q, int v) const
{
  const int z = std::max(_xbar[v], _ybar[v]);
  if (q == _nL)
  {
    return z == 1;
  }
  return _constraints._label[q] == -1 || _constraints._label[q] == z;
}

const IntVector& StateTreeEnumerator::candidates(int q)
{
  IntVector& cand = _cand[q];
  if (_candState[q] == 2)
  {
    return cand;
  }
  cand.clear();
  
  if (_nodeOfIndex[q] != -1)
  {
    cand.push_back(_nodeOfIndex[q]);
    _candState[q] = 2;
    return cand;
  }
  
  int c = -1;
  if (q < _nL)
  {
    c = cell(_states[q].first, _states[q].second);
  }
  else if (_constraints._mutationState >= 0)
  {
    const IntPair& xy = _states[_constraints._mutationState];
    c = cell(xy.first, xy.second);
  }
  
  const int p = q == _nL ? _constraints._mutationState : _constraints._parent[q];
  if (c == -1 || p == -1 || _candState[q] == 1)
  {
    // unconstrained (or cyclic constraints): every node of the copy number state(s) with an
    // allowed label
    for (int v = 0; v < static_cast<int>(_x.size()); ++v)
    {
      if ((c == -1 ? isL(v) : _cell[v] == c) && labelAllowed(q, v))
      {
        cand.push_back(v);
      }
    }
    _candState[q] = 2;
    return cand;
  }
  _candState[q] = 1;
  
  if (q == _nL)
  {
    // targets of the mutation arcs of the candidates of the copy number state of the mutation
    const IntVector& pre = candidates(p);
    for (int v : pre)
    {
      for (int a : _outArcs[v])
      {
        if (_arcs[a]._type == MUTATION)
        {
          cand.push_back(_arcs[a]._target);
        }
      }
    }
    _candState[q] = 2;
    return cand;
  }
  
  // nodes of the copy number state of q with an allowed label, reachable through unused copy
  // number states not in L from the parent and the nodes not in L whose nearest ancestor in L
  // is the parent (if the parent is in T), or from the candidates of the parent (if not)
  const IntVector* parentCandidates = NULL;
  if (p >= 0 && _nodeOfIndex[p] == -1)
  {
    parentCandidates = &candidates(p);
  }
  
  ++_stampNode;
  _queueNode.clear();
  if (parentCandidates == NULL)
  {
    const int pNode = p == -2 ? -1 : _nodeOfIndex[p];
    for (int v : _tNodes)
    {
      if (_owner[v] == pNode)
      {
        _seenNode[v] = _stampNode;
        _queueNode.push_back(v);
      }
    }
  }
  else
  {
    for (int v : *parentCandidates)
    {
      _seenNode[v] = _stampNode;
      _queueNode.push_back(v);
    }
  }
  for (size_t h = 0; h < _queueNode.size(); ++h)
  {
    for (int uw : _outArcs[_queueNode[h]])
    {
      if (_arcs[uw]._type == MUTATION) continue;
      const int w = _arcs[uw]._target;
      const int cw = _cell[w];
      if (_cellUsed[cw] != 0 || _seenNode[w] == _stampNode) continue;
      _seenNode[w] = _stampNode;
      if (_lIndex[cw] < 0)
      {
        _queueNode.push_back(w);
      }
      else if (cw == c && labelAllowed(q, w))
      {
        cand.push_back(w);
      }
    }
  }
  
  _candState[q] = 2;
  return cand;
}

bool StateTreeEnumerator::consistent(int uv) const
{
  if (!_witnessMode) return true;

  const int u = _arcs[uv]._source;
  const int v = _arcs[uv]._target;

  if (_arcs[uv]._type == MUTATION)
  {
    return _constraints._mutationState < 0 || _constraints._mutationState == _lIndex[_cell[v]];
  }

  const int i = _lIndex[_cell[v]];
  if (i < 0) return true;

  const int ownerIndex = _owner[u] == -1 ? -2 : nodeIndex(_owner[u]);
  if (_constraints._parent[i] != -1 && _constraints._parent[i] != ownerIndex)
  {
    return false;
  }
  if (_constraints._label[i] != -1 && _constraints._label[i] != std::max(_xbar[v], _ybar[v]))
  {
    return false;
  }
  return true;
}

void StateTreeEnumerator::report()
{
  ++_stats._reports;

  // short-circuit nodes not in L (StateGraph::finalize)
  StateGraph::StateEdgeSet S;
  for (int v : _tNodes)
  {
    if (!isL(v) || _parent[v] == -1) continue;
    const int p = _owner[_parent[v]];
    if (p == -1) continue;
    S.insert(std::make_pair(triple(p), triple(v)));
  }

  bool ok = false;
  for (const StateGraph::StateEdge& st : S)
  {
    ok |= st.second._z > 0;
  }
  if (!ok) return;

  if (_witnessMode)
  {
    _found = true;
    _witness = S;
    _witnessParents = IntVector(_nL + 1, -2);
    _witnessLabels = IntVector(_nL + 1, 0);
    for (int v : _tNodes)
    {
      if (!isL(v)) continue;
      const int i = nodeIndex(v);
      if (_parent[v] != -1 && _owner[_parent[v]] != -1)
      {
        _witnessParents[i] = nodeIndex(_owner[_parent[v]]);
      }
      _witnessLabels[i] = std::max(_xbar[v], _ybar[v]);
    }
    _witnessMutationState = _lIndex[_cell[_arcs[_mutationArc]._source]];
  }
  else
  {
    _result->insert(S);
  }
}

int StateTreeEnumerator::targetDistance(int uv) const
{
  const int t = _arcs[uv]._target;
  if (_lIndex[_cell[t]] >= 0)
  {
    return consistent(uv) ? 0 : std::numeric_limits<int>::max();
  }
  
  // uncovered states of L that may become state tree children of the owner of the arc's source,
  // with their required label (or any label)
  const int u = _arcs[uv]._source;
  const int ownerIndex = _owner[u] == -1 ? -2 : nodeIndex(_owner[u]);
  int d = std::numeric_limits<int>::max();
  for (int i = 0; i < _nL; ++i)
  {
    if (_cellUsed[cell(_states[i].first, _states[i].second)] != 0
        || (_constraints._parent[i] != -1 && _constraints._parent[i] != ownerIndex)) continue;
    
    const int maxLabel = std::max(_states[i].first, _states[i].second);
    for (int z = 0; z <= maxLabel; ++z)
    {
      if (_constraints._label[i] == -1 || _constraints._label[i] == z)
      {
        d = std::min(d, _distToTarget[_targetOffset[i] + z][t]);
      }
    }
  }
  return d;
}

void StateTreeEnumerator::grow(int depth)
{
  ++_stats._calls;
  
  if (_growBudget != 0 && ++_growCalls > _growBudget)
  {
    _undecided = true;
  }
  if (interrupted())
  {
    return;
  }

  if (_coveredL == _nL && _mutationArc != -1)
  {
    // every further node would be a leaf not in L
    if (_nonLLeaves == 0)
    {
      report();
    }
    return;
  }

  if (_nonLLeaves > _nL - _coveredL)
  {
    ++_stats._prunedDeficit;
    return;
  }

  IntVector& F = _F[depth];

  if (_mutationArc == -1 && !mutationAvailable(F))
  {
    ++_stats._prunedMutation;
    return;
  }

  if (!reachable())
  {
    ++_stats._prunedReach;
    return;
  }

  IntVector& newF = _F[depth + 1];
  IntVector& FF = _FF[depth];
  FF.clear();
  while (!F.empty() && !_found && !interrupted())
  {
    if (_witnessMode)
    {
      // any frontier arc may be picked next; for a fast witness, pick the one whose target is
      // closest to an uncovered state of L
      size_t best = F.size() - 1;
      int bestDist = targetDistance(F[best]);
      for (size_t k = F.size() - 1; k-- > 0 && bestDist > 0;)
      {
        const int d = targetDistance(F[k]);
        if (d < bestDist)
        {
          best = k;
          bestDist = d;
        }
      }
      std::swap(F[best], F.back());
    }
    
    // pick last arc (u,v) from the frontier
    const int uv = F.back();
    F.pop_back();

    if (consistent(uv))
    {
      const int v = _arcs[uv]._target;
      const int cv = _cell[v];

      add(uv);

      // remove arcs into v or into the copy number state of v, and other mutation arcs if (u,v)
      // is the mutation arc
      newF.clear();
      for (int st : F)
      {
        const int t = _arcs[st]._target;
        if (t == v || _cell[t] == cv || (_mutationArc == uv && _cellUsed[_cell[t]] == 1))
          continue;
        newF.push_back(st);
      }

      // push each arc (v,w) where w not in T and the copy number state of w is unused
      // (or (v,w) is a mutation arc and there is no mutation yet)
      for (int vw : _outArcs[v])
      {
        if (_disabled[vw]) continue;
        const int w = _arcs[vw]._target;
        if (_inT[w]) continue;
        const int cw = _cell[w];
        if ((_mutationArc == -1 && cw == cv) || _cellUsed[cw] == 0)
        {
          newF.push_back(vw);
        }
      }

      grow(depth + 1);

      remove(uv);
    }
    else
    {
      ++_stats._prunedConstraint;
    }

    _disabled[uv] = 1;
    FF.push_back(uv);
  }

  for (IntVector::const_reverse_iterator it = FF.rbegin(); it != FF.rend(); ++it)
  {
    F.push_back(*it);
    _disabled[*it] = 0;
  }
}

void StateTreeEnumerator::enumerate(StateGraph::StateEdgeSetSet& result)
{
  _witnessMode = false;
  _undecided = false;
  _growBudget = 0;
  _result = &result;
  init();
  grow(0);
  _result = NULL;
}

bool StateTreeEnumerator::shortestPath(int i, IntVector& path, std::mt19937* rng)
{
  const int q = _constraints._parent[i];
  const int c = cell(_states[i].first, _states[i].second);
  
  ++_stampNode;
  _queueNode.clear();
  const int qNode = q < 0 ? -1 : _nodeOfIndex[q];
  for (int v : _tNodes)
  {
    if (q == -1 || _owner[v] == qNode)
    {
      _seenNode[v] = _stampNode;
      _predArc[v] = -1;
      // state tree parent of a node of L reached from v
      _bfsOwnerIndex[v] = _owner[v] == -1 ? -2 : nodeIndex(_owner[v]);
      _queueNode.push_back(v);
    }
  }
  
  int found = -1;
  for (size_t h = 0; h < _queueNode.size() && found == -1; ++h)
  {
    const int u = _queueNode[h];
    const IntVector& out = _outArcs[u];
    // random tie-breaking between shortest paths: start at a random out-arc
    const size_t offset = rng == NULL || out.empty() ? 0 : (*rng)() % out.size();
    for (size_t k = 0; k < out.size(); ++k)
    {
      const int uw = out[(k + offset) % out.size()];
      if (_arcs[uw]._type == MUTATION || _disabled[uw]) continue;
      const int w = _arcs[uw]._target;
      const int cw = _cell[w];
      if (_cellUsed[cw] != 0 || _seenNode[w] == _stampNode) continue;
      _seenNode[w] = _stampNode;
      _predArc[w] = uw;
      _bfsOwnerIndex[w] = _bfsOwnerIndex[u];
      if (_lIndex[cw] < 0)
      {
        _queueNode.push_back(w);
      }
      else if (cw == c && (q == -1 || q == _bfsOwnerIndex[u]) && labelAllowed(i, w))
      {
        found = w;
        break;
      }
    }
  }
  if (found == -1)
  {
    return false;
  }
  
  // arcs from T to the node found; the path must not repeat a copy number state
  path.clear();
  ++_stamp;
  for (int w = found; _predArc[w] != -1; w = _arcs[_predArc[w]]._source)
  {
    if (_seen[_cell[w]] == _stamp)
    {
      return false;
    }
    _seen[_cell[w]] = _stamp;
    path.push_back(_predArc[w]);
  }
  std::reverse(path.begin(), path.end());
  return true;
}

bool StateTreeEnumerator::route(int allele, int strategy, unsigned int seed)
{
  init();
  std::mt19937 rng(seed);
  std::mt19937* randomness = strategy == 2 ? &rng : NULL;
  
  const int m = _constraints._mutationState;
  if (m < 0)
  {
    return false;
  }
  
  IntVector path;
  IntVector bestPath;
  while (_coveredL < _nL || _mutationArc == -1)
  {
    if (_mutationArc == -1 && _nodeOfIndex[m] != -1)
    {
      // mutate the copy number state of the mutation as soon as it is placed
      int mutation = -1;
      for (int a : _outArcs[_nodeOfIndex[m]])
      {
        if (_arcs[a]._type == MUTATION && (allele == 0 ? _xbar[_arcs[a]._target] : _ybar[_arcs[a]._target]) == 1)
        {
          mutation = a;
        }
      }
      if (mutation == -1 || !consistent(mutation))
      {
        return false;
      }
      add(mutation);
      continue;
    }
    
    // attach an unplaced state of L whose parent is placed, by a shortest path from the
    // region of its parent: the one with the shortest (strategy 0) or longest (1) such path, or
    // a random one (2)
    int best = -1;
    int nReady = 0;
    for (int i = 0; i < _nL; ++i)
    {
      const int q = _constraints._parent[i];
      if (_nodeOfIndex[i] != -1 || (q >= 0 && _nodeOfIndex[q] == -1)) continue;
      if (!shortestPath(i, path, randomness)) continue;
      ++nReady;
      if (best == -1
          || (strategy == 0 && path.size() < bestPath.size())
          || (strategy == 1 && path.size() > bestPath.size())
          || (strategy == 2 && rng() % nReady == 0))
      {
        best = i;
        bestPath = path;
      }
    }
    if (best == -1)
    {
      return false;
    }
    for (int uv : bestPath)
    {
      add(uv);
    }
  }
  
  if (_nonLLeaves != 0)
  {
    return false;
  }
  report();
  return _found;
}

bool StateTreeEnumerator::negotiate(int allele)
{
  const int m = _constraints._mutationState;
  if (m < 0)
  {
    return false;
  }
  for (int i = 0; i < _nL; ++i)
  {
    if (_constraints._parent[i] == -1)
    {
      return false;
    }
  }
  
  init();
  
  // states of L in an order in which every parent precedes its children (the mutated copy, index
  // nL, directly follows m)
  IntVector order;
  std::vector<char> ordered(_nL + 1, 0);
  const int rootIndex = isL(_root) ? _lIndex[_cell[_root]] : -1;
  bool progress = true;
  while (progress)
  {
    progress = false;
    for (int i = 0; i <= _nL; ++i)
    {
      if (ordered[i]) continue;
      const int q = i == _nL ? m : _constraints._parent[i];
      if (i == rootIndex || q == -2 || (q >= 0 && ordered[q]))
      {
        if (i == rootIndex && q != -2) return false;
        ordered[i] = 1;
        order.push_back(i);
        progress = true;
      }
    }
  }
  if (static_cast<int>(order.size()) != _nL + 1)
  {
    return false;
  }
  
  // the arc of the mutation from the (unmutated) node of m, given that node
  std::function<int(int)> mutationArcFrom = [&](int v) {
    for (int a : _outArcs[v])
    {
      if (_arcs[a]._type == MUTATION && (allele == 0 ? _xbar[_arcs[a]._target] : _ybar[_arcs[a]._target]) == 1)
        return a;
    }
    return -1;
  };
  
  // negotiated congestion: route every state of L by a cheapest path from its parent's node,
  // where a copy number state costs more the more other paths use it and the more often it was
  // congested before, until no copy number state is used twice
  IntVector node(_nL + 1, -1);
  IntMatrix pathArcs(_nL + 1);
  IntVector occupancy(_nCells, 0);
  std::vector<double> history(_nCells, 0);
  std::vector<double> dist(_x.size());
  double presentFactor = 0.5;
  typedef std::pair<double, int> QueueItem;
  
  for (int iteration = 0; iteration < NEGOTIATION_ITERATIONS; ++iteration)
  {
    for (int i : order)
    {
      for (int a : pathArcs[i])
      {
        if (!isL(_arcs[a]._target)) --occupancy[_cell[_arcs[a]._target]];
      }
      pathArcs[i].clear();
      
      if (i == rootIndex)
      {
        node[i] = _root;
        continue;
      }
      const int q = i == _nL ? m : _constraints._parent[i];
      const int start = q == -2 ? _root : node[q];
      if (i == _nL)
      {
        const int a = mutationArcFrom(start);
        if (a == -1) return false;
        pathArcs[i].push_back(a);
        node[i] = _arcs[a]._target;
        continue;
      }
      
      // Dijkstra over nodes, through copy number states not in L, to a node of i with an
      // allowed label
      const int c = cell(_states[i].first, _states[i].second);
      std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::infinity());
      std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem> > queue;
      dist[start] = 0;
      _predArc[start] = -1;
      queue.push(QueueItem(0, start));
      int target = -1;
      while (!queue.empty())
      {
        const QueueItem item = queue.top();
        queue.pop();
        const int u = item.second;
        if (item.first > dist[u]) continue;
        if (_cell[u] == c)
        {
          target = u;
          break;
        }
        for (int uw : _outArcs[u])
        {
          if (_arcs[uw]._type == MUTATION) continue;
          const int w = _arcs[uw]._target;
          const int cw = _cell[w];
          double cost;
          if (cw == c)
          {
            if (!labelAllowed(i, w)) continue;
            cost = 1;
          }
          else
          {
            if (_lIndex[cw] >= 0 || cw == _cell[_root]) continue;
            cost = (1 + history[cw]) * (1 + presentFactor * occupancy[cw]);
          }
          if (dist[u] + cost < dist[w])
          {
            dist[w] = dist[u] + cost;
            _predArc[w] = uw;
            queue.push(QueueItem(dist[w], w));
          }
        }
      }
      if (target == -1)
      {
        return false;
      }
      for (int w = target; w != start; w = _arcs[_predArc[w]]._source)
      {
        pathArcs[i].push_back(_predArc[w]);
      }
      std::reverse(pathArcs[i].begin(), pathArcs[i].end());
      node[i] = target;
      for (int a : pathArcs[i])
      {
        if (!isL(_arcs[a]._target)) ++occupancy[_cell[_arcs[a]._target]];
      }
    }
    
    bool congested = false;
    for (int c = 0; c < _nCells; ++c)
    {
      if (occupancy[c] > 1)
      {
        congested = true;
        history[c] += occupancy[c] - 1;
      }
    }
    if (!congested)
    {
      // build the tree with the search's bookkeeping, checking every attachment
      for (int i : order)
      {
        for (int a : pathArcs[i])
        {
          const int w = _arcs[a]._target;
          if (_inT[w] || !_inT[_arcs[a]._source] || (_cellUsed[_cell[w]] != 0 && _arcs[a]._type != MUTATION)
              || !consistent(a))
          {
            return false;
          }
          add(a);
        }
      }
      if (_coveredL != _nL || _mutationArc == -1 || _nonLLeaves != 0)
      {
        return false;
      }
      report();
      return _found;
    }
    presentFactor *= 1.3;
  }
  
  return false;
}

bool StateTreeEnumerator::findWitness(const Constraints& constraints,
                                      StateGraph::StateEdgeSet& witness,
                                      IntVector& parents,
                                      IntVector& labels,
                                      int& mutationState)
{
  _witnessMode = true;
  _constraints = constraints;
  _undecided = false;
  _growCalls = 0;
  
  // the relaxations of the search at the empty tree refute most infeasible constraints at once
  init();
  if (_nonLLeaves > _nL - _coveredL || !mutationAvailable(_F[0]) || !reachable())
  {
    _witnessMode = false;
    return false;
  }
  
  // fast incomplete constructions first (both alleles; shortest-first, longest-first, random
  // orders with random tie-breaking, then negotiated congestion), then the exhaustive search
  bool routed = false;
  for (int attempt = 0; attempt < 2 + ROUTE_RANDOM_ATTEMPTS && !routed; ++attempt)
  {
    const int strategy = std::min(attempt, 2);
    routed = route(0, strategy, attempt) || route(1, strategy, attempt);
  }
  for (int allele = 0; allele < 2 && !routed; ++allele)
  {
    routed = negotiate(allele);
  }
  if (routed)
  {
    ++_stats._routed;
  }
  else
  {
    init();
    grow(0);
  }
  _witnessMode = false;

  if (_found)
  {
    witness = _witness;
    parents = _witnessParents;
    labels = _witnessLabels;
    mutationState = _witnessMutationState;
  }
  return _found;
}
