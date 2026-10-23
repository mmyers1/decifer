/*
 * generatestatetreesmain.cpp
 *
 *  Created on: 28-mar-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "stategraph.h"
#include "readmatrix.h"
#include "arg_parser.h"
#include <fstream>
#include <chrono>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>
#include "solver.h"
#include "enumerator.h"
#include "flashlight.h"
#include "frontier.h"

std::vector<IntPairSet> parseCopyStateInput(std::istream& in)
{
  std::vector<IntPairSet> res;
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of(";"));
    
    IntPairSet L;
    for (const std::string& str : s)
    {
      StringVector ss;
      boost::split(ss, str, boost::is_any_of(","));
      
      if (ss.size() != 2)
      {
        throw std::runtime_error(getLineNumber() + "Error: Invalid copy-number state '" + str + "'");
      }
      
      try
      {
        int x = boost::lexical_cast<int>(ss[0]);
        int y = boost::lexical_cast<int>(ss[1]);
      
      L.insert(IntPair(x, y));
      }
      catch (boost::exception& e)
      {
        throw std::runtime_error(getLineNumber() + "Error: Invalid copy-number state '" + str + "'");
      }
    }
    res.push_back(L);
  }
  
  return res;
}

/// Compute the state trees of L with the 'tree', 'flashlight', 'frontier' or 'auto' algorithm
/// (thread safe)
///
/// @param L Copy number states
/// @param algorithm Algorithm
/// @param frontierLimit State limit per layer of the frontier dynamic program in 'flashlight' and
/// 'auto' (0: none)
/// @param statistics Enumeration statistics and time
StateGraph::StateEdgeSetSet enumerateStateTrees(const IntPairSet& L,
                                                const std::string& algorithm,
                                                int frontierLimit,
                                                std::string& statistics)
{
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  
  StateGraph::StateEdgeSetSet trees;
  StateTreeEnumerator::Stats s;
  FlashlightEnumerator::Stats fs;
  FrontierEnumerator::Stats ds;
  std::string used = algorithm;
  if (algorithm == "auto")
  {
    // run the exact algorithms concurrently and keep the first to finish; their results are
    // identical, so the output does not depend on which one wins. The dynamic program gives up
    // if it exceeds its state limit
    std::atomic<bool> stop(false);
    std::atomic<int> winner(-1);
    StateTreeEnumerator T(L);
    FlashlightEnumerator F(L, frontierLimit);
    FrontierEnumerator D(L);
    T.setStop(&stop);
    F.setStop(&stop);
    D.setStop(&stop);
    D.setStateLimit(frontierLimit);
    StateGraph::StateEdgeSetSet treeTrees;
    StateGraph::StateEdgeSetSet flashlightTrees;
    StateGraph::StateEdgeSetSet frontierTrees;
    std::thread treeThread([&]() {
      T.enumerate(treeTrees);
      if (!stop.exchange(true)) winner = 0;
    });
    std::thread frontierThread([&]() {
      D.enumerate(frontierTrees);
      if (!D.exceeded() && !stop.exchange(true)) winner = 2;
    });
    F.enumerate(flashlightTrees);
    if (!stop.exchange(true)) winner = 1;
    treeThread.join();
    frontierThread.join();
    if (winner == 0)
    {
      trees.swap(treeTrees);
      s = T.stats();
      used = "auto:tree";
    }
    else if (winner == 2)
    {
      trees.swap(frontierTrees);
      ds = D.stats();
      used = "auto:frontier";
    }
    else
    {
      trees.swap(flashlightTrees);
      s = F.searchStats();
      fs = F.stats();
      ds = F.frontierStats();
      used = "auto:flashlight";
    }
  }
  else if (algorithm == "frontier")
  {
    FrontierEnumerator E(L);
    E.enumerate(trees);
    ds = E.stats();
  }
  else if (algorithm == "flashlight")
  {
    FlashlightEnumerator E(L, frontierLimit);
    E.enumerate(trees);
    s = E.searchStats();
    fs = E.stats();
    ds = E.frontierStats();
  }
  else
  {
    StateTreeEnumerator E(L);
    E.enumerate(trees);
    s = E.stats();
  }
  
  std::chrono::duration<double> t = std::chrono::steady_clock::now() - start;
  std::ostringstream out;
  out << " [" << used << " trees=" << trees.size();
  if (used == "flashlight" || used == "auto:flashlight")
  {
    out << " oracle_calls=" << fs._oracleCalls
        << " feasible=" << fs._feasible
        << " reused=" << fs._reused
        << " cycles=" << fs._cycles
        << " routed=" << s._routed
        << " pruned_mandatory=" << s._prunedMandatory
        << " frontier_calls=" << fs._frontierCalls
        << " frontier_exceeded=" << fs._frontierExceeded
        << " frontier_states=" << ds._states
        << " frontier_max_layer=" << ds._maxLayer;
  }
  if (used == "frontier" || used == "auto:frontier")
  {
    out << " states=" << ds._states
        << " partials=" << ds._partials
        << " max_layer=" << ds._maxLayer
        << " max_partials=" << ds._maxPartials
        << " accepted=" << ds._accepted
        << " time=" << t.count() << "s]";
    statistics = out.str();
    return trees;
  }
  out << " calls=" << s._calls
      << " reports=" << s._reports
      << " pruned_deficit=" << s._prunedDeficit
      << " pruned_reach=" << s._prunedReach
      << " pruned_mutation=" << s._prunedMutation
      << " pruned_constraint=" << s._prunedConstraint
      << " time=" << t.count() << "s]";
  statistics = out.str();
  
  return trees;
}

/// Return L with the maternal and paternal copy numbers swapped
IntPairSet swapped(const IntPairSet& L)
{
  IntPairSet res;
  for (const IntPair& xy : L)
  {
    res.insert(IntPair(xy.second, xy.first));
  }
  return res;
}

/// Return the state trees with the maternal and paternal copy numbers swapped. The state graph
/// is symmetric under this swap, so these are the state trees of the swapped copy number states.
StateGraph::StateEdgeSetSet swapped(const StateGraph::StateEdgeSetSet& trees)
{
  StateGraph::StateEdgeSetSet res;
  for (const StateGraph::StateEdgeSet& S : trees)
  {
    StateGraph::StateEdgeSet SS;
    for (const StateGraph::StateEdge& st : S)
    {
      SS.insert(StateGraph::StateEdge(StateGraph::CnaTriple(st.first._y, st.first._x, st.first._z),
                                      StateGraph::CnaTriple(st.second._y, st.second._x, st.second._z)));
    }
    res.insert(SS);
  }
  return res;
}

/// Compute the state trees of L (unless already present) with the given algorithm
///
/// @param L Copy number states
/// @param maxCopy Maximum copy number (ignored by StateGraph::getStateTrees)
/// @param algorithm Algorithm
/// @param stats Whether to report statistics
/// @param symmetry Whether to derive the state trees from those of the swapped copy number
/// states if these are present
/// @param frontierLimit See enumerateStateTrees
void computeStateTrees(const IntPairSet& L,
                       int maxCopy,
                       const std::string& algorithm,
                       bool stats,
                       bool symmetry,
                       int frontierLimit)
{
  if (StateGraph::hasStateTrees(L))
  {
    return;
  }
  
  if (symmetry && algorithm != "legacy" && StateGraph::hasStateTrees(swapped(L)))
  {
    StateGraph::addStateTrees(L, swapped(StateGraph::getStateTrees(swapped(L), maxCopy)));
    if (stats)
    {
      std::cerr << " [derived from the swapped copy number states]";
    }
    return;
  }
  
  if (algorithm == "legacy")
  {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    StateGraph::getStateTrees(L, maxCopy);
    if (stats)
    {
      std::chrono::duration<double> t = std::chrono::steady_clock::now() - start;
      std::cerr << " [legacy trees=" << StateGraph::getStateTrees(L, maxCopy).size()
                << " time=" << t.count() << "s]";
    }
    return;
  }
  
  std::string statistics;
  StateGraph::addStateTrees(L, enumerateStateTrees(L, algorithm, frontierLimit, statistics));
  if (stats)
  {
    std::cerr << statistics;
  }
}

/// Return the progress label of L, e.g. "(1,1) (2,0)"
std::string label(const IntPairSet& L)
{
  std::ostringstream out;
  bool first = true;
  for (const IntPair& xy : L)
  {
    if (!first)
      out << " ";
    else
      first = false;
    out << "(" << xy.first << "," << xy.second << ")";
  }
  return out.str();
}

/// Compute the state trees of every set of copy number states, sequentially (printing progress
/// as before) or with several threads (printing progress as sets finish). The state trees are
/// stored in the ordered dictionary of StateGraph, so the output does not depend on the order.
///
/// @param Ls Sets of copy number states
/// @param labels Progress label of each set
/// @param maxCopies Maximum copy number of each set
/// @param algorithm Algorithm
/// @param stats Whether to report statistics
/// @param threads Number of threads
/// @param symmetry Whether to derive the state trees of sets whose swapped copy number states are
/// present or computed
/// @param frontierLimit See enumerateStateTrees
void computeAllStateTrees(const std::vector<IntPairSet>& Ls,
                          const StringVector& labels,
                          const IntVector& maxCopies,
                          const std::string& algorithm,
                          bool stats,
                          int threads,
                          bool symmetry,
                          int frontierLimit)
{
  if (threads <= 1 || algorithm == "legacy")
  {
    for (size_t i = 0; i < Ls.size(); ++i)
    {
      std::cerr << labels[i] << "..." << std::flush;
      computeStateTrees(Ls[i], maxCopies[i], algorithm, stats, symmetry, frontierLimit);
      std::cerr << " Done." << std::endl;
    }
    return;
  }
  
  // sets to compute, and sets to derive from the swapped copy number states afterwards
  IntVector todo;
  IntVector derived;
  std::set<IntPairSet> seen;
  for (size_t i = 0; i < Ls.size(); ++i)
  {
    if (StateGraph::hasStateTrees(Ls[i]) || !seen.insert(Ls[i]).second) continue;
    const IntPairSet swappedL = swapped(Ls[i]);
    if (symmetry && swappedL != Ls[i]
        && (StateGraph::hasStateTrees(swappedL) || seen.count(swappedL) > 0))
    {
      derived.push_back(i);
    }
    else
    {
      todo.push_back(i);
    }
  }
  
  std::vector<StateGraph::StateEdgeSetSet> results(Ls.size());
  std::atomic<size_t> next(0);
  std::mutex progressMutex;
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t)
  {
    pool.push_back(std::thread([&]() {
      for (size_t k = next++; k < todo.size(); k = next++)
      {
        const int i = todo[k];
        std::string statistics;
        results[i] = enumerateStateTrees(Ls[i], algorithm, frontierLimit, statistics);
        std::lock_guard<std::mutex> lock(progressMutex);
        std::cerr << labels[i] << "..." << (stats ? statistics : "") << " Done." << std::endl;
      }
    }));
  }
  for (std::thread& thread : pool)
  {
    thread.join();
  }
  
  for (int i : todo)
  {
    StateGraph::addStateTrees(Ls[i], results[i]);
  }
  for (int i : derived)
  {
    // the swapped set was present or computed above
    assert(StateGraph::hasStateTrees(swapped(Ls[i])));
    StateGraph::addStateTrees(Ls[i], swapped(StateGraph::getStateTrees(swapped(Ls[i]), maxCopies[i])));
    std::cerr << labels[i] << "..." << (stats ? " [derived from the swapped copy number states]" : "")
              << " Done." << std::endl;
  }
}

bool next(const int maxXY,
          const int maxCN,
          const int nrCnStates,
          const IntVector& minEnumState,
          const IntVector& maxEnumState,
          IntVector& enumState)
{
  assert(enumState.size() == maxCN);
  
  IntVector newEnumState = enumState;
  
  int i = maxCN - 1;
  for (; i >= 0; --i)
  {
    if (newEnumState[i] < maxEnumState[i] - 1)
    {
      ++newEnumState[i];
      break;
    }
  }
  
  if (i < 0)
  {
    return false;
  }
  else
  {
    for (int j = i + 1; j < maxCN; ++j)
    {
      newEnumState[j] = newEnumState[j-1] + 1;
    }
  }
  
  enumState = newEnumState;
  
  return true;
}

int main(int argc, char** argv)
{
  int maxXY = 2;
  int maxCN = 2;
  std::string inputStateTreeFilename;
  std::string outputStateTreeFilename;
  std::string algorithm = "auto";
  bool stats = false;
  int threads = 1;
  bool symmetry = false;
  int frontierLimit = FlashlightEnumerator::DEFAULT_FRONTIER_LIMIT;
  
  lemon2::ArgParser ap(argc, argv);
  ap.refOption("maxXY", "Maximum number of maternal/paternal copies (default: 2)", maxXY, false)
    .refOption("maxCN", "Maximum number of copy number events (default: 2)", maxCN, false)
    .refOption("S", "Input state tree file", inputStateTreeFilename, false)
    .refOption("SS", "Output state tree file", outputStateTreeFilename, false)
    .refOption("algorithm", "Enumeration algorithm: 'tree' (pruned enumeration of state graph subtrees), 'flashlight' (enumeration of state trees with a feasibility oracle), 'frontier' (dynamic program over the lattice of copy number states), 'auto' (tree, flashlight and frontier concurrently, first to finish; uses three threads per set) or 'legacy' (original enumeration); all give identical output (default: auto)", algorithm, false)
    .refOption("stats", "Report enumeration statistics and time per set of copy number states", stats, false)
    .refOption("threads", "Number of threads; sets of copy number states are processed in parallel (not with algorithm 'legacy') (default: 1)", threads, false)
    .refOption("symmetry", "Derive the state trees of a set of copy number states from those of the set with maternal and paternal copy numbers swapped, if that set is also requested (not with algorithm 'legacy')", symmetry, false)
    .refOption("frontierLimit", "Maximum number of states per step of the frontier dynamic program in algorithms 'flashlight' (which then uses its exhaustive search instead) and 'auto' (which then waits for the others); bounds its memory use, about 0.5 kB per state; 0: no limit (default: 2000000)", frontierLimit, false);
  ap.parse();
  
  if (algorithm != "auto" && algorithm != "tree" && algorithm != "flashlight" && algorithm != "frontier" && algorithm != "legacy")
  {
    std::cerr << "Error: invalid algorithm '" << algorithm << "'" << std::endl;
    return 1;
  }
  
  if (!inputStateTreeFilename.empty())
  {
    std::ifstream inS(inputStateTreeFilename.c_str());
    if (!inS.good())
    {
      std::cerr << "Error: could not open '" << inputStateTreeFilename << "' for reading." << std::endl;
      return 1;
    }
    else
    {
      try
      {
        StateGraph::readStateTrees(inS);
      }
      catch (std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
    }
  }
  
  if (!ap.files().empty())
  {
    for (const std::string& filename : ap.files())
    {
      std::ifstream in(filename.c_str());
      if (!in.good())
      {
        std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
        return 1;
      }
      
      std::vector<IntPairSet> cnStates;
      try
      {
        cnStates = parseCopyStateInput(in);
      }
      catch (std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
      in.close();
      
      StringVector labels;
      IntVector maxCopies;
      for (const IntPairSet& L : cnStates)
      {
        int maxCopy = 0;
        for (const IntPair& xy : L)
        {
          maxCopy = std::max(maxCopy, std::max(xy.first, xy.second));
        }
        labels.push_back(label(L));
        maxCopies.push_back(maxCopy);
      }
      computeAllStateTrees(cnStates, labels, maxCopies, algorithm, stats, threads, symmetry, frontierLimit);
      
      std::cerr << "Processed " << filename << "..." << std::endl;
    }
  }
  else
  {
    std::vector<IntPair> cnStates;
    for (int x = 0; x <= maxXY; ++x)
    {
      for (int y = x; y <= maxXY; ++y)
      {
        cnStates.push_back(IntPair(x, y));
        if (x != y)
        {
          cnStates.push_back(IntPair(y, x));
        }
      }
    }
    
    IntVector enumState(maxCN);
    for (int i = 0; i < maxCN; ++i)
    {
      enumState[i] = i;
    }
    
    IntVector minEnumState = enumState;
    
    IntVector maxEnumState(maxCN);
    for (int i = 0; i < maxCN; ++i)
    {
      maxEnumState[i] = cnStates.size() - (maxCN - i - 1);
    }
    
    std::vector<IntPairSet> Ls;
    StringVector labels;
    IntVector maxCopies;
    do
    {
      IntPairSet L;
      int maxCopy = 0;
      std::ostringstream out;
      for (int i = 0; i < maxCN; ++i)
      {
        L.insert(cnStates[enumState[i]]);
        maxCopy = std::max(maxCopy,
                           std::max(cnStates[enumState[i]].first,
                                    cnStates[enumState[i]].second));
        
        if (i > 0)
          out << " ";
        const IntPair& xy = cnStates[enumState[i]];
        out << "(" << xy.first << "," << xy.second << ")";
      }
      
      Ls.push_back(L);
      labels.push_back(out.str());
      maxCopies.push_back(maxCopy + 1);
    } while (next(maxXY, maxCN, cnStates.size(), minEnumState, maxEnumState, enumState));
    
    computeAllStateTrees(Ls, labels, maxCopies, algorithm, stats, threads, symmetry, frontierLimit);
  }
  
  if (outputStateTreeFilename.empty())
  {
    StateGraph::writeStateTrees(std::cout);
  }
  else
  {
    std::ofstream out(outputStateTreeFilename.c_str());
    StateGraph::writeStateTrees(out);
    out.close();
  }
  
  return 0;
}
