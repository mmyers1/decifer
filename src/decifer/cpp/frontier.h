/*
 * frontier.h
 *
 *  Enumeration of state trees by a frontier-based (transfer-matrix) dynamic program over the
 *  lattice of copy number states. Produces the same set of state trees as StateGraph::enumerate,
 *  in time that grows with the width of the lattice instead of with the number of paths.
 */

#ifndef FRONTIER_H
#define FRONTIER_H

#include "stategraph.h"
#include "enumerator.h"
#include <cstdint>
#include <atomic>
#include <string>
#include <set>
#include <vector>
#include <unordered_map>

/// Enumerates the state trees of a set L of copy number states.
///
/// A valid tree T of the state graph (see StateTreeEnumerator) uses every copy number state
/// (x, y) at most once, except the copy number state of the mutation, which holds the
/// unmutated node and the mutated node. Seen on the lattice of copy number states, T is a choice
/// of used cells, of the node (label) in every used cell, and of arcs between adjacent used cells.
/// The cells are processed one at a time in row order. The state of the dynamic program
/// describes the processed cells that still have unprocessed neighbors (the frontier): whether
/// they are used, their label, whether they have a parent and a child, their component in the
/// processed part of T and their owner (nearest ancestor in L, if already known), plus, per
/// component, the nodes of L whose state tree parent will be the owner of the component's
/// future parent. Two partial trees with the same state have the same completions, so every
/// state keeps the set of partial state trees (parent and label of each processed node of L) of
/// the partial trees that reach it, instead of the partial trees.
///
/// Every valid tree corresponds to exactly one sequence of decisions (per cell: unused, its
/// label, or the mutation; per pair of adjacent cells: no arc or an arc in either direction), and
/// every accepted sequence is a valid tree, so the result is exact. Each allele of the mutation is
/// a separate run. States are dropped only by necessary conditions of valid trees: a node with a
/// label other than 0 must be reachable from a mutated node; a node not in L without a child
/// needs a path, through cells not in L, to a node of L that may be a state tree child of its
/// owner; and these open ends may not outnumber the nodes of L that can still become their
/// first descendants in L.
///
/// With constraints (mutation state, state tree parents and labels), the decisions follow them.
/// findWitness also keeps a single partial state tree per state and stops at the first valid
/// tree.
class FrontierEnumerator
{
public:
  struct Stats
  {
    Stats()
      : _states(0)
      , _partials(0)
      , _maxLayer(0)
      , _maxPartials(0)
      , _accepted(0)
    {
    }

    /// Dynamic program states, summed over the cells and both alleles
    uint64_t _states;
    /// Partial state trees kept by the states, summed
    uint64_t _partials;
    /// Largest number of states after processing a cell
    uint64_t _maxLayer;
    /// Largest number of partial state trees of a state
    uint64_t _maxPartials;
    /// Accepted (complete) partial state trees, with repetitions
    uint64_t _accepted;
  };

  /// Constructor
  ///
  /// @param L Copy number states that must be present
  explicit FrontierEnumerator(const IntPairSet& L);

  typedef StateTreeEnumerator::Constraints Constraints;

  /// Insert every state tree (that satisfies the constraints, if given) into result
  ///
  /// @param result Set of state trees
  /// @param constraints Constraints (as in StateTreeEnumerator::findWitness), or NULL
  void enumerate(StateGraph::StateEdgeSetSet& result, const Constraints* constraints = NULL);

  /// Search for one valid tree whose state tree satisfies the constraints (exact: returns false
  /// only if there is none, unless stopped)
  ///
  /// @param constraints Constraints (as in StateTreeEnumerator::findWitness)
  /// @param witness State tree of the valid tree found
  /// @param parents State tree parent of every node of the valid tree found (same indexing)
  /// @param labels Label of every node of the valid tree found (same indexing)
  /// @param mutationState Index in L of the mutation copy number state of the tree found
  /// @return Whether a valid tree was found
  bool findWitness(const Constraints& constraints,
                   StateGraph::StateEdgeSet& witness,
                   IntVector& parents,
                   IntVector& labels,
                   int& mutationState);

  const Stats& stats() const
  {
    return _stats;
  }

  /// Abort the enumeration (with an incomplete result) once *stop is true
  void setStop(const std::atomic<bool>* stop)
  {
    _stop = stop;
  }

  /// Return whether the enumeration has been asked to stop
  bool stopped() const
  {
    return _stop != NULL && _stop->load(std::memory_order_relaxed);
  }

  /// Give up (with an incomplete result, and exceeded() true) when a layer of the dynamic
  /// program has more than this many states (0: no limit). Bounds the memory use
  void setStateLimit(uint64_t limit)
  {
    _stateLimit = limit;
  }

  /// Whether the last enumerate or findWitness gave up because of the state limit
  bool exceeded() const
  {
    return _exceeded;
  }

private:
  /// Frontier cell of a state
  struct Slot
  {
    /// USED, MUT, HAS_PARENT, HAS_CHILD
    uint8_t _flags;
    /// Label (copies of the mutated allele carrying the mutation); the unmutated and mutated
    /// nodes of the mutation cell have labels 0 and 1
    uint8_t _label;
    /// Index in L of the nearest ancestor in L (the node itself if in L), or OWNER_UNKNOWN, or
    /// OWNER_NONE; unused for the mutation cell
    uint8_t _owner;
    /// Component
    uint8_t _comp;
  };

  /// Decoded state, being extended by the decisions of one cell
  struct Work
  {
    bool _mutPlaced;
    /// Component containing the root (-1: root not processed)
    int _rootedComp;
    /// Frontier cells before the current cell, then the current cell
    std::vector<Slot> _slots;
    /// Per component: nodes of L (bit q) whose state tree parent is the owner of the parent of
    /// the component's top
    std::vector<uint64_t> _pending;
    /// Changes to the partial state trees: (index in the partial state tree, value)
    std::vector<std::pair<int, uint8_t> > _delta;
  };

  /// Neighbor of the current cell on the frontier
  struct StepNeighbor
  {
    /// Index in the frontier before the current cell
    int _slot;
    /// Index of the current cell among the neighbors of this cell
    int _toCurrent;
    /// Index of this cell among the neighbors of the current cell
    int _fromCurrent;
  };

  /// Layer of the dynamic program: states (numbered in order of creation) and their partial
  /// state trees
  struct Layer
  {
    std::unordered_map<std::string, int> _index;
    /// Key of each state (in _index)
    std::vector<const std::string*> _keys;
    /// Partial state trees of state e: _partials[_offset[e]], ..., _partials[_offset[e + 1] - 1]
    std::vector<std::string> _partials;
    IntVector _offset;
  };

  int cell(int x, int y) const
  {
    return x * (_maxY + 1) + y;
  }

  /// Return whether some non-mutation arc of the state graph goes from the node of cell c with
  /// label l to the node of its j-th neighbor with label m
  bool arc(int c, int j, int l, int m) const
  {
    return (_arcMask[c][j][l] >> m) & 1;
  }

  void run(int allele, std::set<std::string>& complete);

  /// Whether the search must return early (stop request or state limit)
  bool interrupted() const
  {
    return stopped() || _exceeded;
  }

  /// Processing order of the cells, and the frontiers and step tables that follow from it
  void initOrder();

  void decode(const std::string& key, int k, Work& w) const;

  /// Decide the arcs between the current cell and its frontier neighbors e, e + 1, ...
  void decideArcs(int k, size_t e, const Work& w);

  /// Make the node in slot child a child of the node in slot parent (owner o of the parent's port);
  /// returns false if this gives a node of L a state tree parent that the constraints forbid
  bool attach(Work& w, int child, int parent, uint8_t o) const;

  /// Whether the state tree parent of node q may be p (index, or PARENT_NONE)
  bool parentAllowed(int q, uint8_t p) const;

  /// Whether node q may have label l
  bool labelAllowed(int q, int l) const
  {
    return !_constrained || _constraints._label[q] == -1 || _constraints._label[q] == l;
  }

  /// Whether the mutation may occur in the state of L with index i
  bool mutationAllowed(int i) const
  {
    return !_constrained || _constraints._mutationState == -1 || _constraints._mutationState == i;
  }

  /// Compute the arc labels and _labelOK for the current allele (and constraints)
  void initAllele(int allele);

  /// Partial state tree to state tree (and, if given, parents, labels and mutation state in the
  /// indexing of Constraints)
  StateGraph::StateEdgeSet stateTree(const std::string& p, IntVector* parents, IntVector* labels, int* mutationState) const;

  /// Remove the cells leaving the frontier, and store the resulting state (or accept it)
  void finish(int k, const Work& w);

private:
  const IntPairSet _L;
  std::vector<IntPair> _states;
  int _nL;
  int _maxX;
  int _maxY;
  int _nCells;
  int _rootCell;

  /// Index in L of each cell (-1 if not in L)
  IntVector _lIndex;
  /// Neighbors of each cell that can reach L (empty for cells that cannot)
  IntMatrix _nbrs;
  /// Cells that can reach L
  std::vector<char> _kept;
  /// Cells that can reach L, in processing order
  IntVector _order;
  /// Position of each cell in _order (-1 if absent)
  IntVector _pos;
  /// Largest position of a cell of L
  int _lastL;
  /// _frontier[k]: cells processed before step k with a neighbor processed at step k or later,
  /// in processing order
  IntMatrix _frontier;
  /// Frontier neighbors of the cell of each step
  std::vector<std::vector<StepNeighbor> > _stepNeighbors;
  /// For each step, the index in the work slots (frontier before the step, then the cell of the
  /// step) of every cell of the frontier after the step
  IntMatrix _stepKeep;
  /// For each step, the work slots of the cells that leave the frontier
  IntMatrix _stepLeave;
  /// For each step, whether a cell of L where the mutation may occur is processed later
  std::vector<std::vector<char> > _mutationLater;
  /// For each step, the number of cells of L processed later
  IntVector _lLater;

  /// Current allele (0: maternal, 1: paternal) and the arc labels of the state graph for it:
  /// bit m of _arcMask[c][j][l] (see arc)
  int _allele;
  std::vector<std::vector<std::vector<uint64_t> > > _arcMask;
  /// Whether a node with label l may be in cell c (index c * (max(maxX, maxY) + 1) + l)
  std::vector<char> _labelOK;
  /// _toChild[o][node]: whether a node not in L whose owner is o (index in L, the mutated node
  /// _nL, none _nL + 1, unknown _nL + 2) has a path through cells not in L to a node of L that
  /// may be a state tree child of o, with an allowed label (same node index as _labelOK)
  std::vector<std::vector<char> > _toChild;
  /// Bit o of _allowedParents[q]: whether o may be the state tree parent of q (o = _nL + 1: none)
  std::vector<uint64_t> _allowedParents;

  /// Whether the decisions follow _constraints
  bool _constrained;
  Constraints _constraints;
  /// Existence only (findWitness): one partial state tree per state, stop at the first tree
  bool _existence;
  bool _found;

  // current step: states of the current layer (and the one being processed), the next layer
  // (its partial state trees, by state, are collected in _flat and grouped at the end of the step)
  const Layer* _current;
  size_t _entry;
  Layer _next;
  std::vector<std::pair<int, std::string> > _flat;
  std::set<std::string>* _complete;
  /// Work states by recursion depth (reused, to avoid allocations): 0 decoded, 1 with the
  /// current cell's decision, e + 2 with the arc decision for neighbor e
  std::vector<Work> _pool;
  // scratch space of finish
  std::string _key;
  std::vector<uint8_t> _newComp;
  std::vector<uint64_t> _newPending;

  Stats _stats;
  const std::atomic<bool>* _stop;
  uint64_t _stateLimit;
  bool _exceeded;
};

#endif // FRONTIER_H
