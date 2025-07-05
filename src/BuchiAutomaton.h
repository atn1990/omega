#ifndef OMEGA_BUCHI_AUTOMATON_H
#define OMEGA_BUCHI_AUTOMATON_H

#include <cstdint>
#include <memory>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Automaton.h"
#include "RabinAutomaton.h"

namespace omega {

struct DeterminizeOpts {
  bool SwapUpdateCreate{false};
};

class BuchiAutomaton : public Automaton {
public:
  boost::dynamic_bitset<> initial_states;
  boost::dynamic_bitset<> final_states;

  BuchiAutomaton(size_type alphabet = 0, size_type vertices = 0);

  void Init(const TransitionMap& map);
  void Resize();

  bool Empty();

  void Clean();
  void Reachable();
  void ProjectLabel();

  void Minimize();
  void FindCycle(
      graph_t::vertex_descriptor& s,
      std::vector<int>& component,
      std::vector<int_type>& path,
      char* fmt_s,
      char* fmt_v);
  void FindPath(
      std::vector<graph_t::vertex_descriptor>& pred,
      graph_t::vertex_descriptor& v,
      std::vector<int_type>& path,
      char* fmt_s,
      char* fmt_v);

  std::unique_ptr<TransitionMap> Map() const;

  std::unique_ptr<RabinAutomaton> Determinize(
      const TransitionMap& map,
      DeterminizeOpts opts = DeterminizeOpts()) const;
  std::unique_ptr<BuchiAutomaton> RabinScott(const TransitionMap& map) const;

  friend std::unique_ptr<BuchiAutomaton> Intersection(
      const BuchiAutomaton& A, const BuchiAutomaton& B);
  friend std::unique_ptr<BuchiAutomaton> DisjointUnion(
      const BuchiAutomaton& A, const BuchiAutomaton& B);

  void Print() const;
};

} // namespace omega

#endif // OMEGA_BUCHI_AUTOMATON_H
