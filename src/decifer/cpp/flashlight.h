/*
 * flashlight.h
 *
 *  Enumeration of state trees by backtracking over the decisions that define a state tree,
 *  with a feasibility oracle (flashlight search). Produces the same set of state trees as
 *  StateGraph::enumerate without enumerating the (exponentially many) trees of the state graph
 *  that reduce to the same state tree.
 */

#ifndef FLASHLIGHT_H
#define FLASHLIGHT_H

#include "enumerator.h"
#include "frontier.h"

/// Enumerates the state trees of a set L of copy number states.
///
/// A state tree is defined by the decisions
///  1. the copy number state m of L in which the mutation occurs,
///  2. the parent of every copy number state of L (the mutated copy of m is a child of m),
///  3. the label z of every descendant of the mutated copy (all other labels are 0).
/// Each state tree has exactly one decision sequence. Decisions are made in this order, and a
/// partial decision sequence is extended only if the oracle StateTreeEnumerator::findWitness
/// finds a valid tree consistent with it (or the witness of the previous call already is).
/// The oracle is sound (it returns a valid tree) and complete, so the result equals that of
/// StateTreeEnumerator::enumerate. Its exhaustive search (StateTreeEnumerator::findWitness) is
/// limited to GROW_BUDGET search nodes. If that does not decide, the frontier dynamic program
/// enumerates every state tree consistent with the decisions (FrontierEnumerator::enumerate
/// with constraints), in time that does not grow with the number of paths, and the decisions
/// below are not made: for a parent decision, every state tree with that parent; for a label
/// decision, every state tree with the decided parents (all labels at once). If the dynamic
/// program exceeds its state limit (memory), the decisions are made with the unlimited
/// exhaustive search instead.
class FlashlightEnumerator
{
public:
  struct Stats
  {
    Stats()
      : _oracleCalls(0)
      , _feasible(0)
      , _reused(0)
      , _cycles(0)
      , _frontierCalls(0)
      , _frontierExceeded(0)
    {
    }

    uint64_t _oracleCalls;
    uint64_t _feasible;
    uint64_t _reused;
    uint64_t _cycles;
    /// Oracle calls answered by the frontier dynamic program (subtrees of decisions solved)
    uint64_t _frontierCalls;
    /// Dynamic programs that exceeded the state limit
    uint64_t _frontierExceeded;
  };

  /// Constructor
  ///
  /// @param L Copy number states that must be present
  /// @param frontierLimit State limit per layer of the frontier dynamic program (0: none)
  explicit FlashlightEnumerator(const IntPairSet& L, uint64_t frontierLimit = DEFAULT_FRONTIER_LIMIT);

  /// Insert every state tree into result
  ///
  /// @param result Set of state trees
  void enumerate(StateGraph::StateEdgeSetSet& result);

  const Stats& stats() const
  {
    return _stats;
  }

  /// Abort the enumeration (with an incomplete result) once *stop is true
  void setStop(const std::atomic<bool>* stop)
  {
    _E.setStop(stop);
    _D.setStop(stop);
  }

  /// Statistics of the oracle searches (summed)
  const StateTreeEnumerator::Stats& searchStats() const
  {
    return _E.stats();
  }

  /// Statistics of the dynamic programs of the oracle (summed)
  const FrontierEnumerator::Stats& frontierStats() const
  {
    return _D.stats();
  }

  /// Search nodes of the oracle's exhaustive search before the frontier dynamic program decides
  static const uint64_t GROW_BUDGET = 10000;

  /// Default state limit per layer of the frontier dynamic program (about 1 GB of memory at most)
  static const uint64_t DEFAULT_FRONTIER_LIMIT = 2000000;

private:
  struct Witness
  {
    StateGraph::StateEdgeSet _S;
    IntVector _parents;
    IntVector _labels;
    int _mutationState;
  };

  typedef StateTreeEnumerator::Constraints Constraints;

  typedef enum { INFEASIBLE, FEASIBLE, UNDECIDED } OracleResult;

  /// Return whether the constraints have a valid tree (FEASIBLE, with a valid tree in W), have
  /// none (INFEASIBLE), or whether the budget of the search did not suffice (UNDECIDED; never
  /// if unbounded)
  OracleResult oracle(const Constraints& C, Witness& W, bool unbounded = false);

  /// Insert the state trees of all valid trees that satisfy the constraints into the result, by
  /// the frontier dynamic program; returns false (inserting nothing) if it exceeds its state limit
  bool solve(const Constraints& C);

  /// Return whether the decided parents contain a cycle through node i
  bool cycle(const Constraints& C, int i) const;

  void decideParent(int i, Constraints& C, const Witness& W);

  /// Returns true (at once) if an oracle call was undecided
  bool decideLabel(int k, Constraints& C, const Witness& W);

  /// Make the decisions below the decision with the given oracle result
  void decided(OracleResult r, int i, Constraints& C, Witness& W);

private:
  StateTreeEnumerator _E;
  FrontierEnumerator _D;
  const int _nL;
  /// Index of the normal copy number state (1,1) in L, or -1
  int _rootIndex;
  /// Nodes whose labels are decided, in order
  IntVector _labelNodes;
  /// Whether the label decisions use the unlimited exhaustive search (after the dynamic program
  /// for all labels exceeded its state limit)
  bool _unbounded;
  StateGraph::StateEdgeSetSet* _result;
  Stats _stats;
};

#endif // FLASHLIGHT_H
