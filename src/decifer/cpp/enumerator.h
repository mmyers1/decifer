/*
 * enumerator.h
 *
 *  Exact state tree enumeration with pruning of partial trees that cannot become valid.
 *  Produces the same set of state trees as StateGraph::enumerate.
 */

#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include "stategraph.h"
#include <cstdint>
#include <atomic>
#include <random>

/// Enumerates the state trees of a set L of copy number states.
///
/// The search is the frontier enumeration of StateGraph::grow over rooted subtrees T of the
/// state graph (nodes (x, y, xbar, ybar), root (1,1,0,0)), in which every copy number state
/// (x, y) occurs at most once except for the copy number state of the mutation. A tree T is
/// valid if it contains one mutation edge whose source copy number state is in L, every copy
/// number state of L, and only leaves whose copy number state is in L. The state tree of T is
/// obtained by short-circuiting the nodes whose copy number state is not in L
/// (StateGraph::finalize).
///
/// Partial trees are pruned when no extension can be valid. Every pruning rule is monotone
/// (once violated, it stays violated in every extension), so the set of state trees is exactly
/// that of StateGraph::enumerate:
///  - mutation edges only from copy number states in L
///  - nodes whose copy number state cannot reach any copy number state of L are removed
///  - once L is covered and the mutation is placed, the tree is reported (if all its leaves
///    are in L) and not extended
///  - the number of leaves not in L may not exceed the number of uncovered states of L
///  - every uncovered state of L, and some uncovered state of L from every leaf not in L,
///    must be reachable through unused copy number states
///  - the mutation edge must remain placeable
///
/// In witness mode (findWitness), the search stops at the first valid tree whose state tree is
/// consistent with the given (partial) parents, labels and mutation copy number state. It is
/// further pruned by necessary conditions of the constraints (node-level reachability of every
/// required state tree edge, candidate nodes, and the mandatory cell relaxation), tries greedy
/// constructions first, and can be limited to a budget of search nodes (setGrowBudget).
class StateTreeEnumerator
{
public:
  struct Stats
  {
    Stats()
      : _calls(0)
      , _reports(0)
      , _prunedDeficit(0)
      , _prunedReach(0)
      , _prunedMutation(0)
      , _prunedConstraint(0)
      , _prunedMandatory(0)
      , _routed(0)
    {
    }

    uint64_t _calls;
    uint64_t _reports;
    uint64_t _prunedDeficit;
    uint64_t _prunedReach;
    uint64_t _prunedMutation;
    uint64_t _prunedConstraint;
    /// Pruned by the mandatory cell relaxation (witness mode)
    uint64_t _prunedMandatory;
    /// Witnesses found by route (without the exhaustive search)
    uint64_t _routed;
  };

  /// Constructor
  ///
  /// @param L Copy number states that must be present
  explicit StateTreeEnumerator(const IntPairSet& L);

  /// Insert the state tree of every valid tree into result
  ///
  /// @param result Set of state trees
  void enumerate(StateGraph::StateEdgeSetSet& result);

  /// Constraints of witness mode, indexed by the copy number states of L in IntPairSet order.
  /// The node with index nL denotes the mutated copy of the mutation copy number state.
  struct Constraints
  {
    /// Index in L of the copy number state of the mutation, or -1 if free
    int _mutationState;
    /// Required state tree parent of each node (-1: free, -2: none), size nL + 1
    IntVector _parent;
    /// Required label z of each node (-1: free), size nL + 1
    IntVector _label;
  };

  /// Search for one valid tree whose state tree satisfies the constraints
  ///
  /// @param constraints Constraints
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

  /// Number of copy number states in L
  int nL() const
  {
    return _nL;
  }

  /// Copy number states of L in IntPairSet order
  const std::vector<IntPair>& states() const
  {
    return _states;
  }

  /// Maximum number of maternal copies in the state graph
  int maxX() const
  {
    return _maxX;
  }

  /// Maximum number of paternal copies in the state graph
  int maxY() const
  {
    return _maxY;
  }

  const Stats& stats() const
  {
    return _stats;
  }

  /// Abort the search (enumerate or findWitness return early, with incomplete results) once
  /// *stop is true
  void setStop(const std::atomic<bool>* stop)
  {
    _stop = stop;
  }

  /// Return whether the search has been asked to stop
  bool stopped() const
  {
    return _stop != NULL && _stop->load(std::memory_order_relaxed);
  }

  /// Limit the exhaustive search of findWitness to this many search nodes (0: no limit); when
  /// the limit is reached, findWitness returns false and undecided() returns true
  void setGrowBudget(uint64_t budget)
  {
    _growBudget = budget;
  }

  /// Whether the last findWitness returned false because its search exceeded the budget
  bool undecided() const
  {
    return _undecided;
  }

private:
  typedef enum { MUTATION, AMPLIFICATION, DELETION } ArcType;

  struct StateArc
  {
    int _source;
    int _target;
    ArcType _type;
  };

  int cell(int x, int y) const
  {
    return x * (_maxY + 1) + y;
  }

  bool isL(int v) const
  {
    return _lIndex[_cell[v]] >= 0;
  }

  StateGraph::CnaTriple triple(int v) const
  {
    return StateGraph::CnaTriple(_x[v], _y[v], std::max(_xbar[v], _ybar[v]));
  }

  /// Index of node v in the constraint indexing (copy number states of L, then the mutated copy)
  int nodeIndex(int v) const
  {
    int i = _lIndex[_cell[v]];
    if (i >= 0 && _mutationArc != -1 && _arcs[_mutationArc]._target == v)
      return _nL;
    return i;
  }

  void addArc(int v, ArcType type, int x, int y, int xbar, int ybar);

  void init();

  void add(int uv);

  void remove(int uv);

  bool mutationAvailable(const IntVector& F) const;

  bool reachable();

  bool reachableWitness();

  /// Breadth-first search over nodes from the queue, through unused copy number states not in
  /// L, for a node in the copy number state of L with index target (-1: any uncovered state)
  /// whose constraints allow it to be a state tree child of the node with index ownerIndex
  bool reachesChild(int target, int ownerIndex);

  /// A path that every valid extension of T must contain, through unused copy number states not
  /// in L: from the region of a state of L with a decided parent to it, or from a leaf not in L
  /// to some state of L. All nodes of the path have the same nearest ancestor in L (the owner).
  struct Requirement
  {
    /// Constraint index of the owner of the path's nodes (-2: none)
    int _owner;
    /// Leaf of T the path starts from (-1: it starts in the owner's region)
    int _leaf;
    /// Index in L of the state the path ends in (-1: any uncovered state that may be a state
    /// tree child of the owner)
    int _target;
    /// Nodes the path may start from
    IntVector _sources;
    /// Copy number states on every such path (that avoids the cells reserved for other paths)
    IntVector _mandatory;
  };

  /// Whether the paths of two requirements must be disjoint in every valid extension of T: paths
  /// of different owners (a copy number state not in L has one node, with one owner), and paths
  /// from different leaves (disjoint subtrees)
  static bool disjoint(const Requirement& r, const Requirement& s)
  {
    return r._owner != s._owner || (r._leaf != -1 && s._leaf != -1 && r._leaf != s._leaf);
  }

  /// Breadth-first search for the path of a requirement through unused copy number states not in
  /// L that are not blocked; if found and cells is given, the copy number states of the path
  /// (without its source) are stored in cells
  bool requirementPath(const Requirement& r, IntVector* cells);

  /// Mandatory cell (dominator) relaxation: the copy number states that every path of a
  /// requirement uses, avoiding those mandatory for requirements whose paths must be disjoint
  /// from it, iterated to a fixed point. Returns false if some requirement has no such path
  bool mandatoryCells();

  bool consistent(int uv) const;

  /// Whether the search must return early (stop request or exhausted budget)
  bool interrupted()
  {
    return stopped() || _undecided;
  }

  /// Return whether node v has the label required for the node with constraint index q
  bool labelAllowed(int q, int v) const;

  /// Nodes that the node with constraint index q may be in a valid extension of T: its node if
  /// it is in T; otherwise the nodes of its copy number state with an allowed label that are
  /// reachable, through unused copy number states not in L, from its required parent's region
  /// in T or from its required parent's candidates (computed once per search node)
  const IntVector& candidates(int q);

  /// Smallest number of arcs from the target of arc uv to a node, with an allowed label, of an
  /// uncovered state of L that may become a state tree child of the owner of the arc's source
  /// (0 if the target is in L and consistent with the constraints)
  int targetDistance(int uv) const;

  void grow(int depth);

  /// Shortest path, through unused copy number states not in L, from the region of the required
  /// parent of the state of L with index i (from T if the parent is free) to a node of i that
  /// is consistent with the constraints. Returns false if there is none or if it repeats a
  /// copy number state. Ties between shortest paths are broken at random if rng is given.
  bool shortestPath(int i, IntVector& path, std::mt19937* rng);

  /// Greedy construction of a valid tree satisfying the constraints: repeatedly attach a state
  /// of L whose parent is attached by a shortest path from its parent's region (strategy 0: the
  /// state with the shortest path, 1: the longest, 2: a random one, with random tie-breaking;
  /// seed for 2), and mutate the given allele (0: maternal, 1: paternal) as soon as the copy
  /// number state of the mutation is attached. Incomplete: returns false if it gets stuck. The
  /// tree is built with add() and reported by report(), so it is valid.
  bool route(int allele, int strategy, unsigned int seed);

  /// Number of random routing attempts (per allele) before the exhaustive search
  static const int ROUTE_RANDOM_ATTEMPTS = 8;

  /// Construction of a valid tree satisfying the constraints (all parents must be decided) by
  /// negotiated congestion routing (as in PathFinder, McMurchie and Ebeling 1995): every state
  /// of L is routed from its parent's node by a cheapest path over nodes through copy number
  /// states not in L; copy number states used by several paths become more expensive (present
  /// and history costs) until no copy number state is used twice. Mutates the given allele.
  /// Incomplete: returns false if it does not converge. The tree is built with add() and
  /// reported by report(), so it is valid.
  bool negotiate(int allele);

  /// Number of rounds of negotiate
  static const int NEGOTIATION_ITERATIONS = 200;

  void report();

private:
  const IntPairSet _L;
  std::vector<IntPair> _states;
  int _nL;
  int _maxX;
  int _maxY;
  int _nCells;

  /// Index in L of each copy number state (-1 if not in L)
  IntVector _lIndex;
  /// Copy number states reachable in one amplification or deletion
  IntMatrix _cellNbrs;

  IntVector _x;
  IntVector _y;
  IntVector _xbar;
  IntVector _ybar;
  IntVector _cell;
  IntVector _toNode;
  /// Nodes of each copy number state
  IntMatrix _cellNodes;
  /// Number of arcs from each node to a node of state i of L with label z
  /// (row _targetOffset[i] + z)
  IntMatrix _distToTarget;
  IntVector _targetOffset;
  int _root;

  std::vector<StateArc> _arcs;
  IntMatrix _outArcs;

  // search state
  std::vector<char> _inT;
  std::vector<char> _disabled;
  IntVector _parent;
  IntVector _owner;
  IntVector _childCount;
  IntVector _cellUsed;
  IntVector _tNodes;
  int _nonLLeaves;
  int _coveredL;
  int _mutationArc;

  // frontier and processed arcs per recursion depth
  IntMatrix _F;
  IntMatrix _FF;

  // breadth-first search scratch space
  IntVector _seen;
  int _stamp;
  IntVector _queue;
  IntVector _seenNode;
  int _stampNode;
  IntVector _queueNode;
  IntVector _predArc;
  IntVector _bfsOwnerIndex;
  /// Copy number states the current requirement path may not use
  std::vector<char> _blocked;
  std::vector<Requirement> _requirements;

  /// Node of T with each constraint index (-1 if absent)
  IntVector _nodeOfIndex;
  /// Candidate nodes of each constraint index, and whether they are computed (2), being
  /// computed (1) or not computed (0) at the current search node
  IntMatrix _cand;
  IntVector _candState;

  // witness mode
  bool _witnessMode;
  bool _found;
  Constraints _constraints;
  StateGraph::StateEdgeSet _witness;
  IntVector _witnessParents;
  IntVector _witnessLabels;
  int _witnessMutationState;

  StateGraph::StateEdgeSetSet* _result;
  Stats _stats;
  const std::atomic<bool>* _stop;
  uint64_t _growBudget;
  uint64_t _growCalls;
  bool _undecided;
};

#endif // ENUMERATOR_H
