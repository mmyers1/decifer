/*
 * frontier.cpp
 *
 *  Enumeration of state trees by a frontier-based dynamic program over the lattice of copy
 *  number states.
 */

#include "frontier.h"
#include <algorithm>

namespace {

// Slot flags
const uint8_t USED = 1;
const uint8_t MUT = 2;
const uint8_t HAS_PARENT = 4;
const uint8_t HAS_CHILD = 8;

// Owners: index in L, or
const uint8_t OWNER_NONE = 254;     // no ancestor in L (below the root, which is not in L)
const uint8_t OWNER_UNKNOWN = 255;  // above the top of the component (not yet known)

// Partial state trees: per node q of the state tree (states of L, then the mutated copy),
// byte 2q is the parent (index, PARENT_NONE or UNSET) and byte 2q + 1 the label (or UNSET)
const uint8_t PARENT_NONE = 254;
const uint8_t UNSET = 255;

const uint8_t NO_COMP = 255;

}

FrontierEnumerator::FrontierEnumerator(const IntPairSet& L)
  : _L(L)
  , _states(L.begin(), L.end())
  , _nL(L.size())
  , _maxX(1)
  , _maxY(1)
  , _nCells(0)
  , _rootCell(-1)
  , _lastL(-1)
  , _allele(0)
  , _constrained(false)
  , _existence(false)
  , _found(false)
  , _current(NULL)
  , _entry(0)
  , _complete(NULL)
  , _stop(NULL)
  , _stateLimit(0)
  , _exceeded(false)
{
  // same state graph dimensions as StateGraph::getStateTrees
  for (const IntPair& xy : L)
  {
    _maxX = std::max(_maxX, xy.first);
    _maxY = std::max(_maxY, xy.second);
  }
  // labels are bits of a 64-bit mask, nodes of the state tree bits of a 64-bit pending set
  assert(_maxX < 64 && _maxY < 64 && _nL < 63);
  _nCells = (_maxX + 1) * (_maxY + 1);

  _lIndex = IntVector(_nCells, -1);
  for (int i = 0; i < _nL; ++i)
  {
    _lIndex[cell(_states[i].first, _states[i].second)] = i;
  }

  // copy number moves, as in StateTreeEnumerator, and the cells from which L is reachable (no
  // node of a valid tree is in another cell, since every leaf is in L)
  IntMatrix moves(_nCells);
  for (int x = 0; x <= _maxX; ++x)
  {
    for (int y = 0; y <= _maxY; ++y)
    {
      IntVector& m = moves[cell(x, y)];
      if (0 < x && x < _maxX) m.push_back(cell(x + 1, y));
      if (0 < y && y < _maxY) m.push_back(cell(x, y + 1));
      if (x > 0) m.push_back(cell(x - 1, y));
      if (y > 0) m.push_back(cell(x, y - 1));
    }
  }
  std::vector<char> reachesL(_nCells, 0);
  for (int c = 0; c < _nCells; ++c)
  {
    reachesL[c] = _lIndex[c] >= 0;
  }
  bool changed = true;
  while (changed)
  {
    changed = false;
    for (int c = 0; c < _nCells; ++c)
    {
      if (reachesL[c]) continue;
      for (int d : moves[c])
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
  _rootCell = cell(1, 1);
  assert(reachesL[_rootCell]);

  // lattice neighbors among those cells
  _nbrs = IntMatrix(_nCells);
  for (int x = 0; x <= _maxX; ++x)
  {
    for (int y = 0; y <= _maxY; ++y)
    {
      const int c = cell(x, y);
      if (!reachesL[c]) continue;
      const int dx[] = {1, 0, -1, 0};
      const int dy[] = {0, 1, 0, -1};
      for (int k = 0; k < 4; ++k)
      {
        const int xx = x + dx[k];
        const int yy = y + dy[k];
        if (0 <= xx && xx <= _maxX && 0 <= yy && yy <= _maxY && reachesL[cell(xx, yy)])
        {
          _nbrs[c].push_back(cell(xx, yy));
        }
      }
    }
  }

  _kept = reachesL;
  initOrder();
}

void FrontierEnumerator::initOrder()
{
  // processing order: rows along the shorter side of the lattice, so the frontier is short,
  // from the first row (near the root) up. Sweeping from the last row down (or from the row of
  // the mutation in both directions, which doubles the frontier) was slower in tests
  const bool rows = _maxX <= _maxY;
  const int nRows = (rows ? _maxY : _maxX) + 1;
  const int nCols = (rows ? _maxX : _maxY) + 1;
  IntVector rowOrder;
  for (int r = 0; r < nRows; ++r)
  {
    rowOrder.push_back(r);
  }
  _order.clear();
  _pos = IntVector(_nCells, -1);
  for (int r : rowOrder)
  {
    for (int i = 0; i < nCols; ++i)
    {
      const int c = rows ? cell(i, r) : cell(r, i);
      if (_kept[c]) _order.push_back(c);
    }
  }
  const int n = _order.size();
  for (int k = 0; k < n; ++k)
  {
    _pos[_order[k]] = k;
  }
  _lastL = -1;
  for (int i = 0; i < _nL; ++i)
  {
    _lastL = std::max(_lastL, _pos[cell(_states[i].first, _states[i].second)]);
  }

  // a cell is on the frontier from its step until the step of its last neighbor
  IntVector last(_nCells, -1);
  for (int c : _order)
  {
    last[c] = _pos[c];
    for (int d : _nbrs[c])
    {
      last[c] = std::max(last[c], _pos[d]);
    }
  }
  _frontier = IntMatrix(n + 1);
  for (int k = 0; k <= n; ++k)
  {
    for (int c : _order)
    {
      if (_pos[c] < k && last[c] >= k)
      {
        _frontier[k].push_back(c);
      }
    }
  }

  _stepNeighbors = std::vector<std::vector<StepNeighbor> >(n);
  _stepKeep = IntMatrix(n);
  _stepLeave = IntMatrix(n);
  for (int k = 0; k < n; ++k)
  {
    const int v = _order[k];
    const IntVector& before = _frontier[k];
    const int nBefore = before.size();
    for (int s = 0; s < nBefore; ++s)
    {
      const int u = before[s];
      const IntVector::const_iterator it = std::find(_nbrs[u].begin(), _nbrs[u].end(), v);
      if (it == _nbrs[u].end()) continue;
      StepNeighbor sn;
      sn._slot = s;
      sn._toCurrent = it - _nbrs[u].begin();
      sn._fromCurrent = std::find(_nbrs[v].begin(), _nbrs[v].end(), u) - _nbrs[v].begin();
      _stepNeighbors[k].push_back(sn);
    }
    for (int c : _frontier[k + 1])
    {
      _stepKeep[k].push_back(c == v ? nBefore : std::find(before.begin(), before.end(), c) - before.begin());
    }
    for (int s = 0; s <= nBefore; ++s)
    {
      if (std::find(_stepKeep[k].begin(), _stepKeep[k].end(), s) == _stepKeep[k].end())
      {
        _stepLeave[k].push_back(s);
      }
    }
  }

  // number of cells of L processed after step k
  _lLater = IntVector(n, 0);
  for (int k = 0; k < n; ++k)
  {
    for (int i = 0; i < _nL; ++i)
    {
      if (_pos[cell(_states[i].first, _states[i].second)] > k)
      {
        ++_lLater[k];
      }
    }
  }
  
  // whether the mutation (of each allele) may still occur in a cell processed after step k
  _mutationLater = std::vector<std::vector<char> >(2, std::vector<char>(n, 0));
  for (int a = 0; a < 2; ++a)
  {
    for (int k = 0; k < n; ++k)
    {
      for (int i = 0; i < _nL; ++i)
      {
        const IntPair& xy = _states[i];
        if (_pos[cell(xy.first, xy.second)] > k && (a == 0 ? xy.first : xy.second) > 0)
        {
          _mutationLater[a][k] = 1;
        }
      }
    }
  }
}

void FrontierEnumerator::decode(const std::string& key, int k, Work& w) const
{
  const int n = _frontier[k].size();
  w._mutPlaced = key[0] != 0;
  w._rootedComp = static_cast<uint8_t>(key[1]) == NO_COMP ? -1 : static_cast<uint8_t>(key[1]);
  w._slots.resize(n);
  for (int s = 0; s < n; ++s)
  {
    Slot& slot = w._slots[s];
    slot._flags = key[2 + 4 * s];
    slot._label = key[3 + 4 * s];
    slot._owner = key[4 + 4 * s];
    slot._comp = key[5 + 4 * s];
  }
  const int nComp = (key.size() - 2 - 4 * n) / 8;
  w._pending.resize(nComp);
  for (int c = 0; c < nComp; ++c)
  {
    uint64_t p = 0;
    for (int b = 0; b < 8; ++b)
    {
      p |= static_cast<uint64_t>(static_cast<uint8_t>(key[2 + 4 * n + 8 * c + b])) << (8 * b);
    }
    w._pending[c] = p;
  }
  w._delta.clear();
}

bool FrontierEnumerator::parentAllowed(int q, uint8_t p) const
{
  if (!_constrained || _constraints._parent[q] == -1) return true;
  return _constraints._parent[q] == (p == PARENT_NONE ? -2 : static_cast<int>(p));
}

bool FrontierEnumerator::attach(Work& w, int child, int parent, uint8_t o) const
{
  const uint8_t A = w._slots[parent]._comp;
  const uint8_t B = w._slots[child]._comp;
  assert(A != B && static_cast<int>(B) != w._rootedComp);
  w._slots[child]._flags |= HAS_PARENT;
  w._slots[parent]._flags |= HAS_CHILD;

  // the nodes of L waiting for the owner of the top's parent get it (or keep waiting, in A)
  const uint64_t pending = w._pending[B];
  w._pending[B] = 0;
  for (int q = 0; q <= _nL; ++q)
  {
    if (!((pending >> q) & 1)) continue;
    if (o == OWNER_UNKNOWN)
    {
      w._pending[A] |= static_cast<uint64_t>(1) << q;
    }
    else
    {
      const uint8_t p = o == OWNER_NONE ? PARENT_NONE : o;
      if (!parentAllowed(q, p)) return false;
      w._delta.push_back(std::make_pair(2 * q, p));
    }
  }

  // merge B into A
  for (Slot& s : w._slots)
  {
    if (!(s._flags & USED) || s._comp != B) continue;
    if (s._owner == OWNER_UNKNOWN)
    {
      s._owner = o;
    }
    s._comp = A;
  }
  return true;
}

void FrontierEnumerator::decideArcs(int k, size_t e, const Work& w)
{
  const std::vector<StepNeighbor>& neighbors = _stepNeighbors[k];
  if (e == neighbors.size())
  {
    finish(k, w);
    return;
  }

  // no arc
  decideArcs(k, e + 1, w);

  const int vs = w._slots.size() - 1;
  const StepNeighbor& sn = neighbors[e];
  const Slot& U = w._slots[sn._slot];
  const Slot& V = w._slots[vs];
  if (!(U._flags & USED) || !(V._flags & USED) || U._comp == V._comp) return;
  const int u = _frontier[k][sn._slot];
  const int v = _order[k];

  // arc from (a node of) u to v; the mutated node of a mutation cell has the unmutated one as
  // its parent, so arcs into a mutation cell go to the unmutated node (label 0)
  if (!(V._flags & HAS_PARENT))
  {
    const int target = (V._flags & MUT) ? 0 : V._label;
    for (int port = 0; port < ((U._flags & MUT) ? 2 : 1); ++port)
    {
      const int label = (U._flags & MUT) ? port : U._label;
      const uint8_t owner = (U._flags & MUT) ? (port == 0 ? _lIndex[u] : _nL) : U._owner;
      if (arc(u, sn._toCurrent, label, target))
      {
        Work& w2 = _pool[e + 2];
        w2 = w;
        if (attach(w2, vs, sn._slot, owner))
        {
          decideArcs(k, e + 1, w2);
        }
      }
    }
  }

  // arc from (a node of) v to u
  if (!(U._flags & HAS_PARENT))
  {
    const int target = (U._flags & MUT) ? 0 : U._label;
    for (int port = 0; port < ((V._flags & MUT) ? 2 : 1); ++port)
    {
      const int label = (V._flags & MUT) ? port : V._label;
      const uint8_t owner = (V._flags & MUT) ? (port == 0 ? _lIndex[v] : _nL) : V._owner;
      if (arc(v, sn._fromCurrent, label, target))
      {
        Work& w2 = _pool[e + 2];
        w2 = w;
        if (attach(w2, sn._slot, vs, owner))
        {
          decideArcs(k, e + 1, w2);
        }
      }
    }
  }
}

void FrontierEnumerator::finish(int k, const Work& w)
{
  const int nBefore = _frontier[k].size();

  // cells leaving the frontier get no further arcs: every node needs a parent (the root has one
  // by convention) and every node not in L a child
  for (int s : _stepLeave[k])
  {
    const Slot& slot = w._slots[s];
    if (!(slot._flags & USED)) continue;
    if (!(slot._flags & HAS_PARENT)) return;
    const int c = s < nBefore ? _frontier[k][s] : _order[k];
    if (!(slot._flags & MUT) && _lIndex[c] < 0 && !(slot._flags & HAS_CHILD)) return;
  }

  // components on the new frontier, numbered in order of appearance
  const IntVector& keep = _stepKeep[k];
  std::vector<uint8_t>& newComp = _newComp;
  newComp.assign(w._pending.size(), NO_COMP);
  int nComp = 0;
  for (int s : keep)
  {
    const Slot& slot = w._slots[s];
    if ((slot._flags & USED) && newComp[slot._comp] == NO_COMP)
    {
      newComp[slot._comp] = nComp++;
    }
  }

  // a component that left the frontier is complete. Only the one containing the root can be
  // (the top of any other one has no parent); the tree is then complete if nothing else is on
  // the frontier, every state of L is processed and the mutation is placed
  for (int s : _stepLeave[k])
  {
    const Slot& slot = w._slots[s];
    if ((slot._flags & USED) && newComp[slot._comp] == NO_COMP)
    {
      if (slot._comp == w._rootedComp && nComp == 0 && w._mutPlaced && k >= _lastL)
      {
        for (int i = _current->_offset[_entry]; i < _current->_offset[_entry + 1]; ++i)
        {
          std::string q = _current->_partials[i];
          for (const std::pair<int, uint8_t>& d : w._delta)
          {
            q[d.first] = d.second;
          }
          ++_stats._accepted;
          _complete->insert(q);
          if (_existence)
          {
            _found = true;
            break;
          }
        }
      }
      return;
    }
  }

  if (!w._mutPlaced && !_mutationLater[_allele][k])
  {
    return;
  }
  
  // a node not in L without a child (an open end) needs a path to an allowed state tree child of
  // its owner. Open ends of a component have disjoint subtrees, so each needs its own first
  // descendant in L: an unprocessed state of L or a pending node of another component. Two open
  // ends can share one only if the component of one is attached below the other, which happens
  // once per component without the root (as R4 in StateTreeEnumerator)
  const int K = std::max(_maxX, _maxY) + 1;
  int openEnds = 0;
  for (size_t i = 0; i < keep.size(); ++i)
  {
    const Slot& slot = w._slots[keep[i]];
    if (!(slot._flags & USED) || (slot._flags & (MUT | HAS_CHILD))) continue;
    const int c = _frontier[k + 1][i];
    if (_lIndex[c] >= 0) continue;
    const int o = slot._owner == OWNER_UNKNOWN ? _nL + 2 : (slot._owner == OWNER_NONE ? _nL + 1 : slot._owner);
    if (!_toChild[o][c * K + slot._label]) return;
    ++openEnds;
  }
  int available = _lLater[k];
  int attachments = 0;
  for (size_t c = 0; c < w._pending.size(); ++c)
  {
    if (newComp[c] == NO_COMP || static_cast<int>(c) == w._rootedComp) continue;
    ++attachments;
    for (int q = 0; q <= _nL; ++q)
    {
      available += (w._pending[c] >> q) & 1;
    }
  }
  if (openEnds - attachments > available)
  {
    return;
  }
  
  // the nodes of L waiting for the same owner need a common allowed state tree parent
  if (_constrained)
  {
    for (size_t c = 0; c < w._pending.size(); ++c)
    {
      if (newComp[c] == NO_COMP || w._pending[c] == 0) continue;
      uint64_t allowed = ~static_cast<uint64_t>(0);
      for (int q = 0; q <= _nL; ++q)
      {
        if ((w._pending[c] >> q) & 1)
        {
          allowed &= _allowedParents[q];
        }
      }
      if (allowed == 0) return;
    }
  }

  std::string& key = _key;
  key.clear();
  key.push_back(w._mutPlaced ? 1 : 0);
  key.push_back(w._rootedComp == -1 ? NO_COMP : newComp[w._rootedComp]);
  for (size_t i = 0; i < keep.size(); ++i)
  {
    const Slot& slot = w._slots[keep[i]];
    if (slot._flags & USED)
    {
      // whether a node has a child matters only for nodes not in L
      const bool inL = (slot._flags & MUT) || _lIndex[_frontier[k + 1][i]] >= 0;
      key.push_back(inL ? slot._flags & ~HAS_CHILD : slot._flags);
      key.push_back(slot._label);
      key.push_back(slot._owner);
      key.push_back(newComp[slot._comp]);
    }
    else
    {
      key.append(4, 0);
    }
  }
  std::vector<uint64_t>& pending = _newPending;
  pending.assign(nComp, 0);
  for (size_t c = 0; c < newComp.size(); ++c)
  {
    if (newComp[c] != NO_COMP)
    {
      pending[newComp[c]] = w._pending[c];
    }
  }
  for (uint64_t p : pending)
  {
    for (int b = 0; b < 8; ++b)
    {
      key.push_back(static_cast<char>((p >> (8 * b)) & 0xFF));
    }
  }

  std::unordered_map<std::string, int>::iterator it = _next._index.find(key);
  if (it == _next._index.end())
  {
    it = _next._index.insert(std::make_pair(key, static_cast<int>(_next._keys.size()))).first;
    _next._keys.push_back(&it->first);
    if (_stateLimit != 0 && _next._keys.size() > _stateLimit)
    {
      _exceeded = true;
    }
  }
  else if (_existence)
  {
    // existence only: one partial tree per state suffices
    return;
  }
  const int e = it->second;
  for (int i = _current->_offset[_entry]; i < _current->_offset[_entry + 1]; ++i)
  {
    _flat.push_back(std::make_pair(e, _current->_partials[i]));
    for (const std::pair<int, uint8_t>& d : w._delta)
    {
      _flat.back().second[d.first] = d.second;
    }
    if (_existence) break;
  }
}

void FrontierEnumerator::initAllele(int allele)
{
  _allele = allele;
  
  // labels of the arcs of the state graph (as in StateTreeEnumerator) between nodes whose label
  // counts copies of this allele (the other allele has no mutated copy)
  _arcMask = std::vector<std::vector<std::vector<uint64_t> > >(_nCells);
  for (int c : _order)
  {
    const int x = c / (_maxY + 1);
    const int y = c % (_maxY + 1);
    const int maxLabel = allele == 0 ? x : y;
    _arcMask[c] = std::vector<std::vector<uint64_t> >(_nbrs[c].size(), std::vector<uint64_t>(maxLabel + 1, 0));
    for (int l = 0; l <= maxLabel; ++l)
    {
      const int xbar = allele == 0 ? l : 0;
      const int ybar = allele == 0 ? 0 : l;
      std::vector<std::vector<int> > targets;
      if (0 < x && x < _maxX && xbar < x) targets.push_back({x + 1, y, xbar, ybar});
      if (0 < x && x < _maxX && xbar > 0) targets.push_back({x + 1, y, xbar + 1, ybar});
      if (0 < y && y < _maxY && ybar < y) targets.push_back({x, y + 1, xbar, ybar});
      if (0 < y && y < _maxY && ybar > 0) targets.push_back({x, y + 1, xbar, ybar + 1});
      if (x > xbar) targets.push_back({x - 1, y, xbar, ybar});
      if (xbar > 0) targets.push_back({x - 1, y, xbar - 1, ybar});
      if (y > ybar) targets.push_back({x, y - 1, xbar, ybar});
      if (ybar > 0) targets.push_back({x, y - 1, xbar, ybar - 1});
      for (const std::vector<int>& t : targets)
      {
        const int d = cell(t[0], t[1]);
        const IntVector::const_iterator it = std::find(_nbrs[c].begin(), _nbrs[c].end(), d);
        if (it == _nbrs[c].end()) continue;
        assert((allele == 0 ? t[3] : t[2]) == 0);
        _arcMask[c][it - _nbrs[c].begin()][l] |= static_cast<uint64_t>(1) << (allele == 0 ? t[2] : t[3]);
      }
    }
  }
  
  // necessary conditions on the node (c, l) of a valid tree, by breadth-first search over the
  // arcs (ignoring that copy number states are used once): if l > 0, it descends from the
  // mutated node (m, 1) of an allowed mutation state m; if c is not in L, it has a path through
  // cells not in L to a node of L that may be a state tree child of its owner, with an allowed
  // label (_toChild, per owner)
  const int K = std::max(_maxX, _maxY) + 1;
  const int nNodes = _nCells * K;
  IntMatrix out(nNodes);
  IntMatrix in(nNodes);
  for (int c : _order)
  {
    for (size_t j = 0; j < _nbrs[c].size(); ++j)
    {
      for (size_t l = 0; l < _arcMask[c][j].size(); ++l)
      {
        for (int m = 0; m < K; ++m)
        {
          if ((_arcMask[c][j][l] >> m) & 1)
          {
            out[c * K + l].push_back(_nbrs[c][j] * K + m);
            in[_nbrs[c][j] * K + m].push_back(c * K + l);
          }
        }
      }
    }
  }
  std::vector<char> fromMutation(nNodes, 0);
  IntVector queue;
  for (int i = 0; i < _nL; ++i)
  {
    const int c = cell(_states[i].first, _states[i].second);
    const int maxLabel = _allele == 0 ? _states[i].first : _states[i].second;
    if (maxLabel > 0 && mutationAllowed(i) && labelAllowed(_nL, 1))
    {
      fromMutation[c * K + 1] = 1;
      queue.push_back(c * K + 1);
    }
  }
  for (size_t h = 0; h < queue.size(); ++h)
  {
    for (int u : out[queue[h]])
    {
      if (!fromMutation[u])
      {
        fromMutation[u] = 1;
        queue.push_back(u);
      }
    }
  }
  
  // owners: index in L (the mutated node included), or _nL + 1 for none
  _toChild = std::vector<std::vector<char> >(_nL + 3, std::vector<char>(nNodes, 0));
  std::vector<char>& toAny = _toChild[_nL + 2];
  for (int o = 0; o <= _nL + 1; ++o)
  {
    std::vector<char>& to = _toChild[o];
    IntVector queueL;
    for (int i = 0; i < _nL; ++i)
    {
      if (i == o || !parentAllowed(i, o == _nL + 1 ? PARENT_NONE : o)) continue;
      const int c = cell(_states[i].first, _states[i].second);
      const int maxLabel = _allele == 0 ? _states[i].first : _states[i].second;
      for (int l = 0; l <= maxLabel; ++l)
      {
        if (labelAllowed(i, l))
        {
          to[c * K + l] = 1;
          queueL.push_back(c * K + l);
        }
      }
    }
    for (size_t h = 0; h < queueL.size(); ++h)
    {
      for (int u : in[queueL[h]])
      {
        if (!to[u] && _lIndex[u / K] < 0)
        {
          to[u] = 1;
          queueL.push_back(u);
        }
      }
    }
    for (int u = 0; u < nNodes; ++u)
    {
      toAny[u] |= to[u];
    }
  }
  
  // state tree parents each node of L may have (bit _nL + 1: none)
  _allowedParents = std::vector<uint64_t>(_nL + 1, 0);
  for (int q = 0; q <= _nL; ++q)
  {
    for (int o = 0; o <= _nL + 1; ++o)
    {
      if (o != q && parentAllowed(q, o == _nL + 1 ? PARENT_NONE : o))
      {
        _allowedParents[q] |= static_cast<uint64_t>(1) << o;
      }
    }
  }
  
  _labelOK = std::vector<char>(nNodes, 0);
  for (int c : _order)
  {
    const int maxLabel = _allele == 0 ? c / (_maxY + 1) : c % (_maxY + 1);
    for (int l = 0; l <= maxLabel; ++l)
    {
      _labelOK[c * K + l] = (l == 0 || fromMutation[c * K + l]) && (_lIndex[c] >= 0 || toAny[c * K + l]);
    }
  }
}

void FrontierEnumerator::run(int allele, std::set<std::string>& complete)
{
  initAllele(allele);
  _complete = &complete;
  const int K = std::max(_maxX, _maxY) + 1;
  
  Layer current;
  const std::string start = std::string(1, 0) + std::string(1, static_cast<char>(NO_COMP));
  current._keys.push_back(&current._index.insert(std::make_pair(start, 0)).first->first);
  current._partials.push_back(std::string(2 * (_nL + 1), static_cast<char>(UNSET)));
  current._offset = {0, 1};
  
  const int n = _order.size();
  _pool.resize(8);
  Work& w = _pool[0];
  for (int k = 0; k < n && !interrupted() && !_found; ++k)
  {
    const int v = _order[k];
    const int x = v / (_maxY + 1);
    const int y = v % (_maxY + 1);
    const int li = _lIndex[v];
    const bool isRoot = v == _rootCell;
    // the mutation must be in this cell (witness mode with a decided mutation state)
    const bool mustMutate = li >= 0 && _constrained && _constraints._mutationState == li;
    
    _next._index.clear();
    _next._keys.clear();
    _flat.clear();
    _current = &current;
    for (_entry = 0; _entry < current._keys.size() && !interrupted() && !_found; ++_entry)
    {
      decode(*current._keys[_entry], k, w);
      const int comp = w._pending.size();
      
      // unused (not for states of L and the root)
      if (li < 0 && !isRoot)
      {
        w._slots.push_back(Slot());
        w._slots.back()._flags = 0;
        w._slots.back()._label = 0;
        w._slots.back()._owner = 0;
        w._slots.back()._comp = 0;
        finish(k, w);
        w._slots.pop_back();
      }
      
      // one node with label l (the root has label 0 and, by convention, a parent)
      const int maxLabel = isRoot ? 0 : (allele == 0 ? x : y);
      for (int l = 0; l <= maxLabel && !mustMutate; ++l)
      {
        if (!_labelOK[v * K + l] || (li >= 0 && !labelAllowed(li, l))) continue;
        if (isRoot && li >= 0 && !parentAllowed(li, PARENT_NONE)) continue;
        Work& w2 = _pool[1];
        w2 = w;
        Slot slot;
        slot._flags = USED | (isRoot ? HAS_PARENT : 0);
        slot._label = l;
        slot._owner = li >= 0 ? li : (isRoot ? OWNER_NONE : OWNER_UNKNOWN);
        slot._comp = comp;
        w2._slots.push_back(slot);
        w2._pending.push_back(li >= 0 && !isRoot ? static_cast<uint64_t>(1) << li : 0);
        if (isRoot)
        {
          w2._rootedComp = comp;
        }
        if (li >= 0)
        {
          w2._delta.push_back(std::make_pair(2 * li + 1, static_cast<uint8_t>(l)));
          if (isRoot)
          {
            w2._delta.push_back(std::make_pair(2 * li, PARENT_NONE));
          }
        }
        decideArcs(k, 0, w2);
      }
      
      // the mutation: the unmutated node (label 0) and its child, the mutated node (label 1)
      if (li >= 0 && !w._mutPlaced && (allele == 0 ? x : y) > 0 && mutationAllowed(li) && labelAllowed(li, 0)
          && labelAllowed(_nL, 1) && parentAllowed(_nL, li) && (!isRoot || parentAllowed(li, PARENT_NONE)))
      {
        Work& w2 = _pool[1];
        w2 = w;
        Slot slot;
        slot._flags = USED | MUT | (isRoot ? HAS_PARENT : 0);
        slot._label = 0;
        slot._owner = li;
        slot._comp = comp;
        w2._slots.push_back(slot);
        w2._pending.push_back(isRoot ? 0 : static_cast<uint64_t>(1) << li);
        w2._mutPlaced = true;
        if (isRoot)
        {
          w2._rootedComp = comp;
          w2._delta.push_back(std::make_pair(2 * li, PARENT_NONE));
        }
        w2._delta.push_back(std::make_pair(2 * li + 1, static_cast<uint8_t>(0)));
        w2._delta.push_back(std::make_pair(2 * _nL, static_cast<uint8_t>(li)));
        w2._delta.push_back(std::make_pair(2 * _nL + 1, static_cast<uint8_t>(1)));
        decideArcs(k, 0, w2);
      }
    }
    
    // partial state trees, grouped by state; those reaching a state through different partial
    // trees are kept once
    std::sort(_flat.begin(), _flat.end());
    _flat.erase(std::unique(_flat.begin(), _flat.end()), _flat.end());
    const int nStates = _next._keys.size();
    _next._offset.assign(nStates + 1, 0);
    for (const std::pair<int, std::string>& f : _flat)
    {
      ++_next._offset[f.first + 1];
    }
    for (int e = 0; e < nStates; ++e)
    {
      _stats._maxPartials = std::max<uint64_t>(_stats._maxPartials, _next._offset[e + 1]);
      _next._offset[e + 1] += _next._offset[e];
    }
    _next._partials.clear();
    for (std::pair<int, std::string>& f : _flat)
    {
      _next._partials.push_back(std::string());
      _next._partials.back().swap(f.second);
    }
    _stats._states += nStates;
    _stats._partials += _flat.size();
    _stats._maxLayer = std::max<uint64_t>(_stats._maxLayer, nStates);
    
    std::swap(current, _next);
  }
  _next = Layer();
  _flat.clear();
  _current = NULL;
  _complete = NULL;
}

StateGraph::StateEdgeSet FrontierEnumerator::stateTree(const std::string& p,
                                                       IntVector* parents,
                                                       IntVector* labels,
                                                       int* mutationState) const
{
  const int m = static_cast<uint8_t>(p[2 * _nL]);
  StateGraph::StateEdgeSet S;
  for (int q = 0; q <= _nL; ++q)
  {
    const uint8_t parent = p[2 * q];
    const uint8_t label = p[2 * q + 1];
    assert(parent != UNSET && label != UNSET);
    if (parents != NULL)
    {
      (*parents)[q] = parent == PARENT_NONE ? -2 : parent;
      (*labels)[q] = label;
    }
    if (parent == PARENT_NONE) continue;
    const IntPair& xy = _states[parent < _nL ? parent : m];
    const IntPair& uv = _states[q < _nL ? q : m];
    S.insert(StateGraph::StateEdge(StateGraph::CnaTriple(xy.first, xy.second, static_cast<uint8_t>(p[2 * parent + 1])),
                                   StateGraph::CnaTriple(uv.first, uv.second, label)));
  }
  if (mutationState != NULL)
  {
    *mutationState = m;
  }
  return S;
}

bool FrontierEnumerator::findWitness(const Constraints& constraints,
                                     StateGraph::StateEdgeSet& witness,
                                     IntVector& parents,
                                     IntVector& labels,
                                     int& mutationState)
{
  _constrained = true;
  _existence = true;
  _constraints = constraints;
  _found = false;
  std::set<std::string> complete;
  _exceeded = false;
  for (int allele = 0; allele < 2 && !_found && !interrupted(); ++allele)
  {
    run(allele, complete);
  }
  _constrained = false;
  _existence = false;
  if (!_found)
  {
    return false;
  }
  _found = false;
  parents = IntVector(_nL + 1, -2);
  labels = IntVector(_nL + 1, 0);
  witness = stateTree(*complete.begin(), &parents, &labels, &mutationState);
  return true;
}

void FrontierEnumerator::enumerate(StateGraph::StateEdgeSetSet& result, const Constraints* constraints)
{
  _constrained = constraints != NULL;
  _existence = false;
  if (_constrained)
  {
    _constraints = *constraints;
  }
  std::set<std::string> complete;
  _exceeded = false;
  for (int allele = 0; allele < 2 && !interrupted(); ++allele)
  {
    run(allele, complete);
  }
  _constrained = false;
  if (interrupted())
  {
    return;
  }
  
  // state tree of each complete partial state tree
  for (const std::string& p : complete)
  {
    result.insert(stateTree(p, NULL, NULL, NULL));
  }
}
