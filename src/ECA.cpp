/**
  * <ECA.cpp>
  *   Functions to construct and test properties of Büchi automata.
 **/

#include "BuchiAutomaton.h"
#include "ECA.h"

#include <cmath>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <print>
#include <queue>

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/strong_components.hpp>

#define DE_BRUIJN_SIZE        4
#define DE_BRUIJN_TRANSITIONS 8
#define DE_BRUIJN_ALPHABET    2

#define ECA_SIZE     10
#define ECA_ALPHABET 4
#define ECA_INITIAL  0

namespace omega {

// Print the rule table of an elementary cellular automaton.
void RuleTable(int_type rule) {
  std::print("# RuleTable({})\n", rule);

  // 8 transitions in total
  for (auto i = 0U; i < DE_BRUIJN_TRANSITIONS; i++) {
    std::print("{}{}{} = {}\n", map_bit(i, 2), map_bit(i, 1), map_bit(i, 0), map_bit(rule, i));
  }

  std::print("\n");
}

struct GlobalMapOpts {
  bool full = true;
  bool negated = false;
};

// Constructs a Büchi automaton that verifies the global map for a given pair of tracks.
// rule is the ECA rule number
// k is the number of tracks
// pair indicates which tracks evolve under the global map
// opts indicates some options for constructing the automaton
std::unique_ptr<BuchiAutomaton> GlobalMap(int_type rule, int_type k, int_pair pair, GlobalMapOpts opts = GlobalMapOpts()) {
  BOOST_ASSERT_MSG(pair.first < k, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.second < k, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.first != pair.second, "Tracks must be different");

  // S = 2^k, N = 10
  auto num_states = ECA_SIZE;
  if (!opts.full) {
    num_states--;
  }

  auto M = std::make_unique<BuchiAutomaton>(std::exp2(k), num_states);

  dbg(OutputType::General, {
    std::print("# GlobalMap(");
    if (opts.negated) {
      std::print("x{} -/> x{})", pair.first, pair.second);
    } else {
      std::print("x{} --> x{})", pair.first, pair.second);
    }
    std::print("  N = {}, k = {}, S = 2^k = {}\n\n", num_states, k, M->num_alphabet);
  });

  auto type = boost::get(boost::vertex_name, M->graph);
  auto label = boost::get(boost::edge_name, M->graph);

  // If negated == true, then full == true.
  BOOST_ASSERT(opts.negated ? opts.full : true);

  // Crash state (omitted if opts.full == false).
  // If opts.negated == true, then it is a final state.
  graph_t::vertex_descriptor crash;
  if (opts.full) {
    crash = boost::vertex(ECA_SIZE-1, M->graph);

    if (opts.negated) {
      M->final_states.set(ECA_SIZE-1);
      type[crash] = NodeType::Final;
    } else {
      type[crash] = NodeType::None;
    }

    // Absorbing state, so add an edge for every symbol to itself.
    for (auto i = 0UL; i < M->num_alphabet; i++) {
      auto [e, added] = boost::add_edge(crash, crash, M->graph);
      label[e] = i;
    }
  }

  // Initial state.
  auto initial = boost::vertex(0, M->graph);
  type[initial] = NodeType::Initial;
  M->initial_states.set(0);

  auto width = binary_digits(M->num_alphabet-1);
  dbg(OutputType::Debug, std::print("# width = {}\n", width));

  // There are 4 "initial" states of the form 0:x2:y2
  for (auto i = 0U; i < M->num_alphabet; i++) {
    // decompose i into two bits x2:y2
    auto x2 = map_bit(i, pair.first);
    auto y2 = map_bit(i, pair.second);

    // The numbering of states corresponds to the binary representation of its
    // component bits, accounting for the initial state.
    auto u = boost::vertex(compose_3(0, x2, y2) + 1, M->graph);

    dbg(OutputType::Debug, std::print("*  {:0{}b}  -->  [0{} {}]\n", i, width, x2, y2));

    auto [e, b] = boost::add_edge(initial, u, M->graph);
    label[e] = i;
  }

  dbg(OutputType::Debug, std::print("\n"));

  // There are 8 possible combinations of x1:x2:y2 in total
  for (auto i = 0U; i < ECA_SIZE-2; i++) {
    auto u = boost::vertex(i+1, M->graph);

    // decompose i into three bits x1:x2:y2
    auto x1 = get_bit(i, 2);
    auto x2 = get_bit(i, 1);
    auto y2 = get_bit(i, 0);

    // 4 possible combinations of x3:y3
    for (auto j = 0U; j < M->num_alphabet; j++) {
      // decompose j into two bits x3:y3
      auto x3 = get_bit(j, pair.first);
      auto y3 = get_bit(j, pair.second);
      auto v = boost::vertex(compose_3(x2, x3, y3) + 1, M->graph);

      dbg(OutputType::Debug, std::print("{}  [{}{} {}]  {:0{}b}", i, x1, x2, y2, j, width));

      // x1:x2 -> x2:x3
      //   :y2      :y3
      // if and only if MAP(rule, x1:x2:x3) == y2
      if (map_bit(rule, compose_3(x1, x2, x3)) == y2) {
        dbg(OutputType::Debug, std::print("  -->  {}  [{}{} {}]\n", compose_3(x2, x3, y3), x2, x3, y3));
        dbg(OutputType::Debug, std::print("  {}{}{} --> {}\n", x1, x2, x3, y2));

        auto [e, added] = boost::add_edge(u, v, M->graph);
        label[e] = j;
      } else {
        dbg(OutputType::Debug, std::print("  -->  X\n"));

        if (opts.full) {
          auto [e, added] = boost::add_edge(u, crash, M->graph);
          label[e] = j;
        }
      }
    }

    dbg(OutputType::Debug, std::print("\n"));

    if (opts.negated) {
      type[u] = NodeType::None;
    } else {
      type[u] = NodeType::Final;
      M->final_states.set(i+1);
    }
  }

  M->num_vertices = boost::num_vertices(M->graph);
  M->num_edges = boost::num_edges(M->graph);

  dbg(OutputType::Debug, std::print("\n"));
  dbg(OutputType::Debug, M->Print());

  if (!opts.full) {
    M->Clean();

    // Since initial states are also final.
    // TODO: remove if changed.
    // M->final_states.reset(0);

    dbg(OutputType::General, M->Print());
  }

  return M;
}

// Construct a Büchi automaton to check if two tracks are equal.
// k is the number of tracks
// pair indicates which tracks should be checked for equality
std::unique_ptr<BuchiAutomaton> Equality(int_type k, int_pair pair) {
  BOOST_ASSERT_MSG(pair.first < k, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.second < k, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.first != pair.second, "Tracks must be different");
  auto M = std::make_unique<BuchiAutomaton>(std::exp2(k), 1);

  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  auto u = boost::vertex(0, M->graph);
  type[u] = NodeType::Both;

  for (auto i = 0U; i < M->num_alphabet; i++) {
    if (map_bit(i, pair.first) == map_bit(i, pair.second)) {
      auto [e, added] = boost::add_edge(u, u, M->graph);
      BOOST_ASSERT(added);

      label[e] = i;
    }
  }

  M->Resize();
  dbg(OutputType::General, M->Print());

  return M;
}

// Modify a Büchi automaton to ensure that two specified tracks are equal.
void Equality(BuchiAutomaton& M, int_pair pair) {
  BOOST_ASSERT_MSG(static_cast<int_type>(std::exp2(pair.first)) < M.num_alphabet, "Track index out of bounds");
  BOOST_ASSERT_MSG(static_cast<int_type>(std::exp2(pair.second)) < M.num_alphabet, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.first != pair.second, "Tracks must be different");

  dbg(OutputType::General, std::print("# x{} == x{}\n", pair.first, pair.second));

  auto label = boost::get(boost::edge_name, M.graph);

  boost::remove_edge_if([&pair, &label](const auto& edge) {
    auto symbol = label[edge];
    return map_bit(symbol, pair.first) != map_bit(symbol, pair.second);
  }, M.graph);

  M.Resize();
  dbg(OutputType::General, M.Print());

  M.Reachable();
  dbg(OutputType::General, M.Print());
}

// Specialized product construction with the inequality automaton.
// The product P is two copies of M, one where the tracks have not differed,
// and one where they have.
std::unique_ptr<BuchiAutomaton> Inequality(BuchiAutomaton& M, int_pair pair) {
  BOOST_ASSERT_MSG(static_cast<int_type>(std::exp2(pair.first)) < M.num_alphabet, "Track index out of bounds");
  BOOST_ASSERT_MSG(static_cast<int_type>(std::exp2(pair.second)) < M.num_alphabet, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.first != pair.second, "Tracks must be different");

  auto P = std::make_unique<BuchiAutomaton>(M.num_alphabet);

  dbg(OutputType::Debug, std::print("# Inequality({}, {})\n", pair.first, pair.second));
  dbg(OutputType::Debug, std::print("# (N = {}, S = {})\n\n", M.num_vertices, M.num_alphabet));
  dbg(OutputType::Debug, std::print("# Initial States\n"));

  auto vertex_width = decimal_digits(M.num_vertices);
  char fmt_vrtx[MAXLINE];
  snprintf(fmt_vrtx, MAXLINE, "(%%%du, %%d)\n", vertex_width);

  std::unordered_map<std::pair<int_type, bool>, graph_t::vertex_descriptor> map;
  std::queue<std::tuple<int_type, bool, graph_t::vertex_descriptor>> queue;

  auto label = boost::get(boost::edge_name, P->graph);
  auto type = boost::get(boost::vertex_name, P->graph);

  auto label_M = boost::get(boost::edge_name, M.graph);
  auto index_M = boost::get(boost::vertex_index, M.graph);

  auto EQ = false;
  auto NE = true;

  for (auto i : dynamic_bitset_iterator(M.initial_states)) {
    // Add a vertex to the graph representing the state (i, EQ).
    auto u = boost::add_vertex(P->graph);

    // Store (i, EQ) -> u in the map.
    map.insert({{i, EQ}, u});

    // Fill in the graph as a BFS from the initial vertices.
    queue.push({i, EQ, u});

    // The state (i, EQ) is an initial state in P.
    type[u] = NodeType::Initial;

    dbg(OutputType::Debug, printf(fmt_vrtx, i, EQ));
  }

  dbg(OutputType::Debug, std::print("\n"));

  char fmt_edge[MAXLINE];
  snprintf(fmt_edge, MAXLINE, "  -->  (%%%du, %%d)", vertex_width);

  auto edge_width = binary_digits(P->num_alphabet-1);
  dbg(OutputType::Debug, std::print("# edge width = {}\n", edge_width));

  while (!queue.empty()) {
    auto [source, source_component, u] = queue.front();
    queue.pop();

    dbg(OutputType::Debug, printf(fmt_vrtx, source, source_component));

    for (auto [e_itr, e_end] = boost::out_edges(boost::vertex(source, M.graph), M.graph); e_itr != e_end; ++e_itr) {
      auto symbol = label_M[*e_itr];

      auto target = index_M[boost::target(*e_itr, M.graph)];
      auto component = source_component;

      // If the tracks are not equal, go into the second component.
      if (get_bit(symbol, pair.first) != get_bit(symbol, pair.second)) {
        component = NE;
      }

      dbg(OutputType::Debug, {
        std::print("  {:0{}b}", symbol, edge_width);
        printf(fmt_edge, target, component);
      });

      auto itr = map.find({target, component});

      // Add a new vertex to the product machine.
      if (itr == map.end()) {
        dbg(OutputType::Debug, std::print("  *"));

        auto v = boost::add_vertex(P->graph);

        if (component == NE && M.final_states[target]) {
          type[v] = NodeType::Final;
        } else {
          type[v] = NodeType::None;
        }

        itr = map.insert({{target, component}, v}).first;

        queue.push({target, component, v});
      }

      // Add a new edge from u to the target vertex.
      auto [e, added] = boost::add_edge(u, itr->second, P->graph);
      label[e] = symbol;

      dbg(OutputType::Debug, std::print("\n"));
    }

    dbg(OutputType::Debug, std::print("\n"));
  }

  P->Resize();

  // Remove all initial states from the set of final states.
  // P->final_states &= ~P->initial_states;

  dbg(OutputType::General, P->Print());

  return P;
}

// Construct a Büchi automaton to check if two tracks are not equal.
std::unique_ptr<BuchiAutomaton> Inequality(int_type k, int_pair pair) {
  auto M = std::make_unique<BuchiAutomaton>(std::exp2(k), 2);

  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  auto u = boost::vertex(0, M->graph);
  auto v = boost::vertex(1, M->graph);

  type[u] = NodeType::Initial;
  type[v] = NodeType::Final;

  for (auto i = 0U; i < M->num_alphabet; i++) {
    if (get_bit(i, pair.first) == get_bit(i, pair.second)) {
      auto [e, added] = boost::add_edge(u, u, M->graph);
      label[e] = i;
    } else {
      auto [e, added] = boost::add_edge(u, v, M->graph);
      label[e] = i;
    }

    auto [e, added] = boost::add_edge(v, v, M->graph);
    label[e] = i;
  }

  M->num_edges = 2 * M->num_alphabet;
  M->initial_states.set(0);
  M->final_states.set(1);

  dbg(OutputType::General, M->Print());

  return M;
}

// Construct a Büchi automaton to check if track is a forbidden configuration.
std::unique_ptr<BuchiAutomaton> Pattern(int_type k, int_type n, const std::string& s) {
  auto start = 0U;

  graph_t::vertex_descriptor u, v, w;
  graph_t::edge_descriptor e;
  bool added;

  auto M = std::make_unique<BuchiAutomaton>(std::exp2(k), s.length() + 2);
  auto label  = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  // u is the final state
  // if a configuration doesn't match the pattern, it will enter the trap state
  // however, if an arbitrary transient is allowed, the initial state is the
  // final state
  if (s[0] == '2') {
    u = boost::vertex(1, M->graph);
    start = 1;
  } else {
    u = boost::vertex(s.length() + 1, M->graph);
  }

  for (auto i = start; i < s.length(); i++) {
    auto v = boost::vertex(i, M->graph);
    auto w = boost::vertex(i+1, M->graph);

    type[v] = NodeType::None;

    if (s[i] == '*') {
      if (i < s.length() - 1) {
        i++;
        v = boost::vertex(i-2, M->graph);
        w = boost::vertex(i+1, M->graph);
      } else {
        continue;
      }
    }

    for (auto j = 0U; j < M->num_alphabet; j++) {
      int_type symbol = s[i % s.length()] - '0';
      if (map_bit(j, n) == symbol) {
        if (i < s.length()-1 && s[i+1] == '*') {
          std::tie(e, added) = boost::add_edge(v, v, M->graph);
        } else {
          std::tie(e, added) = boost::add_edge(v, w, M->graph);
        }
      } else if (i < s.length()-1 && s[i+1] == '*') {
        continue;
      } else {
        std::tie(e, added) = boost::add_edge(v, u, M->graph);
      }

      label[e] = j;
    }
  }

  if (s[s.length() - 1] != '*') {
    v = boost::vertex(s.length(), M->graph);
  } else {
    v = boost::vertex(s.length() - 2, M->graph);
  }

  if (s[start + 1] != '*') {
    w = boost::vertex(start + 1, M->graph);
  } else {
    w = boost::vertex(start + 3, M->graph);
  }

  type[v] = NodeType::None;

  int_type symbol;
  for (auto j = 0UL; j < M->num_alphabet; j++) {
    if (s[start + 1] != '*') {
      symbol = s[start] - '0';
      if (map_bit(j, n) == symbol) {
        auto [e, b] = boost::add_edge(v, w, M->graph);
      } else if (s[s.length() - 1] != '*') {
        auto [e, b] = boost::add_edge(v, u, M->graph);
      }

      label[e] = j;
    } else {
      symbol = s[0] - '0';
      if (map_bit(j, n) == symbol) {
        auto [e, b] = boost::add_edge(v, v, M->graph);
        label[e] = j;
      }

      if (map_bit(j, n) == (int_type) (s[2] - '0')) {
        auto [e, b] = boost::add_edge(v, w, M->graph);
        label[e] = j;
      } else if (map_bit(j, n) != symbol) {
        auto [e, b] = boost::add_edge(v, u, M->graph);
        label[e] = j;
      }
    }

    if (s[0] != '2') {
      auto [e, b] = boost::add_edge(u, u, M->graph);
      label[e] = j;
    }
  }

  type[u] = NodeType::Final;

  v = boost::vertex(start, M->graph);
  type[v] = NodeType::Initial;

  dbg(OutputType::General, M->Print());

  M->Clean();

  if (s[0] != '2') {
    M->final_states.reset(0);
  }

  dbg(OutputType::General, M->Print());

  return M;
}

// Construct a Büchi automaton that checks if a track is a finite configuration (but has at least one).
std::unique_ptr<BuchiAutomaton> Finite(int_type k, int_type n, const std::string& s) {
  auto num_states = 2;
  auto M = std::make_unique<BuchiAutomaton>(std::exp2(k), num_states);
  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  auto u = boost::vertex(0, M->graph);
  type[u] = NodeType::Initial;

  auto v = boost::vertex(1, M->graph);
  type[v] = NodeType::Final;

  // finite occurrence of s[1]
  for (auto j = 0UL; j < M->num_alphabet; j++) {
    if (map_bit(j, n) != (int_type) (s[1] - '0')) {
      auto [e, b] = boost::add_edge(v, v, M->graph);
      label[e] = j;
    }

    auto [e, b] = boost::add_edge(u, u, M->graph);
    label[e] = j;

    std::tie(e, b) = boost::add_edge(u, v, M->graph);
    label[e] = j;
  }

  M->Resize();
  M->final_states.reset(0);

  dbg(OutputType::General, M->Print());

  auto N = std::make_unique<BuchiAutomaton>(std::exp2(k), 2);
  label = boost::get(boost::edge_name, N->graph);
  type = boost::get(boost::vertex_name, N->graph);

  u = boost::vertex(0, N->graph);
  type[u] = NodeType::Initial;

  v = boost::vertex(1, N->graph);
  type[v] = NodeType::Final;

  // at least one occurrence of s[1]
  for (auto j = 0UL; j < N->num_alphabet; j++) {
    if (map_bit(j, n) != (int_type) (s[1] - '0')) {
      auto [e, b] = boost::add_edge(u, u, N->graph);
      label[e] = j;
    }

    else {
      auto [e, b] = boost::add_edge(u, v, N->graph);
      label[e] = j;
    }

    auto [e, b] = boost::add_edge(v, v, N->graph);
    label[e] = j;
  }

  N->Resize();
  N->final_states.reset(0);

  dbg(OutputType::General, N->Print());

  return Intersection(*M, *N);
}

std::string to_string(ShiftType shift) {
  switch (shift) {
    case ShiftType::Left:
      return "Left";
    case ShiftType::Right:
      return "Right";
  }
}

// Construct a Büchi automaton that checks if track1 is a left or right shift of track2.
// @param k: the number of tracks
// @param pair: a pair of indices indicating which tracks to compare
// @param shift: a ShiftType indicating whether to perform a left or right shift
// @param full: if true, the automaton will include a crash state that absorbs all transitions
//              if false, the automaton will not include the crash state
// @return: a unique pointer to the constructed Büchi automaton
// The automaton will have 4 states if `full` is true, and 3 states otherwise.
// The states are:
// 0: initial state
// 1: final state for track1
// 2: final state for track2
// 3: crash state (only if `full` is true)
// The transitions are defined such that:
// - From the initial state, transition to either state 1 or 2 based on the corresponding track's symbol.
// - From state 1, transition to state 1 or 2 if the corresponding track's symbol matches the previous symbol read.
// - If the corresponding track's symbol does not match, transition to the crash state if `full` is true.
// - Similarly, from state 2, transition to state 1 or 2 based on the corresponding track's symbol, or to the crash state if `full` is true and the symbol does not match.
// - The crash state absorbs all transitions, meaning any symbol will loop back to itself.
std::unique_ptr<BuchiAutomaton> Shift(int_type k, int_pair pair, ShiftType shift, bool full = false) {
  BOOST_ASSERT_MSG(pair.first < k, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.second < k, "Track index out of bounds");
  BOOST_ASSERT_MSG(pair.first != pair.second, "Tracks must be different");

  dbg(OutputType::Debug, std::print("# {}-Shift({}, {{{}, {}}})\n", to_string(shift), k, pair.first, pair.second));

  auto num_states = 4;
  if (!full) {
    num_states--;
  }

  auto M = std::make_unique<BuchiAutomaton>(std::exp2(k), num_states);

  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  // initial state
  auto initial = boost::vertex(0, M->graph);
  type[initial] = NodeType::Initial;
  M->initial_states.set(0);

  // crash state
  auto crash = initial;
  if (full) {
    crash = boost::vertex(3, M->graph);
    type[crash] = NodeType::None;
  }

  for (auto i = 0UL; i < M->num_alphabet; i++) {
    auto n = 0;
    if (shift == ShiftType::Left) {
      n = get_bit(i, pair.second);
    } else {
      n = get_bit(i, pair.first);
    }

    auto u = boost::vertex(n+1, M->graph);

    // add edge from initial state
    auto [e, added] = boost::add_edge(initial, u, M->graph);
    label[e] = i;

    if (full) {
      // add edge on crash state
      auto [e, added] = boost::add_edge(crash, crash, M->graph);
      label[e] = i;
    }
  }

  for (auto i = 0UL; i < 2; i++) {
    auto u = boost::vertex(i+1, M->graph);
    type[u] = NodeType::Final;

    for (auto j = 0UL; j < M->num_alphabet; j++) {
      auto n = 0;
      if (shift == ShiftType::Left) {
        n = get_bit(j, pair.second);
      } else {
        n = get_bit(j, pair.first);
      }
      auto v = boost::vertex(n+1, M->graph);

      // previous symbol of x matches the current symbol of y
      if (shift == ShiftType::Left) {
        if (get_bit(j, pair.first) == i) {
          auto [e, added] = boost::add_edge(u, v, M->graph);
          label[e] = j;
        } else if (full) {
          auto [e, added] = boost::add_edge(u, crash, M->graph);
          label[e] = j;
        }
      } else {
        if (get_bit(j, pair.second) == i) {
          auto [e, added] = boost::add_edge(u, v, M->graph);
          label[e] = j;
        } else if (full) {
          auto [e, added] = boost::add_edge(u, crash, M->graph);
          label[e] = j;
        }
      }
    }
  }

  M->Resize();

  dbg(OutputType::General, M->Print());

  return M;
}

bool Shift(int_type rule, int_type k, ShiftType shift) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k >= 1, "Parameter k must be at least 1");

  dbg(OutputType::General, std::print("# Shift({}, {}, {})\n", rule, k, to_string(shift)));

  GlobalMapOpts opts;
  opts.full = false;

  auto M = GlobalMap(rule, k+1, {0, 1}, opts);

  for (auto i = 1U; i < k; i++) {
    auto N = GlobalMap(rule, k+1, {i, i+1}, opts);

    dbg(OutputType::Debug, std::print("# x0 -> x{}\n", i+1));
    M = Intersection(*M, *N);
  }

  for (auto i = 0U; i < k; i++) {
    M = Inequality(*M, {i, i+1});
  }

  auto S = Shift(k+1, {0, k}, shift);
  M = Intersection(*M, *S);

  // for (const auto& s : p) {
  //   dbg(OutputType::General, std::print("# x != {}\n", s));
  //   N = Pattern(k+1, 0, s);

  //   dbg(OutputType::General, std::print("# [x -> y] && x != {}\n", s));
  //   M = Intersection(*M, *N);

  //   std::print("{}  ", !M->Empty());
  // }

  auto result = !M->Empty();
  dbg(OutputType::Quiet, std::print("Shift({}, {}, {}): {}\n", rule, k, to_string(shift), result));
  return result;
}

// Construct a Büchi automaton to check if an ECA has a fixed-point.
// A one-way infinite elementary cellular automaton has a fixed-point if there exists a configuration that evolves to itself after one application of the global map.
bool FixedPoint(int_type rule) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, std::print("# FixedPoint({})\n", rule));

  // Construct x_0 -> x_1 and x_0 == x_1
  dbg(OutputType::General, std::print("# x0 -> x1\n"));
  auto M = GlobalMap(rule, 2, {0, 1}, opts);

  dbg(OutputType::General, std::print("# x0 == x1\n"));
  Equality(*M, {0, 1});

  auto result = !M->Empty();
  dbg(OutputType::Quiet, std::print("FixedPoint({}): {}\n", rule, result));
  return result;

  // std::unique_ptr<BuchiAutomaton> N;
  // ensure fixed point isn't forbidden
  // for (auto i = 0UL; i < p.size(); i++) {
  //   // finite support
  //   if (p[i][0] == 'f') {
  //     if ((p[i][1] == '1') && MAP(rule, 0)) {
  //       std::print("X  ");
  //       continue;
  //     }

  //     else if ((p[i][1] == '0') && !MAP(rule, 7)) {
  //       std::print("X  ");
  //       continue;
  //     }

  //     dbg(OutputType::General, std::print("# x0 == {}\n", p[i]));
  //     N = Finite(2, 0, p[i]);
  //     dbg(OutputType::General, std::print("# [x -> x] && x0 == {}\n", p[i]));
  //   } else {
  //     dbg(OutputType::General, std::print("# x0 != {}\n", p[i]));
  //     N = Pattern(2, 0, p[i]);
  //     dbg(OutputType::General, std::print("# [x -> x] && x0 != {}\n", p[i]));
  //   }

  //   M = Intersection(*M, *N);

  //   // std::print("{}  ", !M->Empty());
  //   return !M->Empty();
  // }
}

// Construct a Büchi automaton to check if an ECA has a k-cycle.
//
// A one-way infinite elementary cellular automaton has a k-cycle if there exists a configuration that evolves to itself after k applications of the global map.
//
// Furthermore, the cycle must be proper, i.e., not a d-cycle where d is a divisor of k.
bool Cycle(int_type rule, int_type k) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k >= 1, "Parameter k must be at least 1");

  GlobalMapOpts opts;
  opts.full = false;

  // construct x_1 -> x_2
  // for i = {2, ..., k-1}, construct x_i -> x_{i+1} and x_1 -> x_{i+1}
  // finally, construct x_k -> x_1 and x_1 -> x_1 from x_1 -> x_k

  dbg(OutputType::General, std::print("# Cycle({}, {})\n", rule, k));
  auto M = GlobalMap(rule, k+1, {0, 1}, opts);

  for (auto i = 1U; i < k; i++) {
    auto N = GlobalMap(rule, k+1, {i, i+1}, opts);

    dbg(OutputType::General, std::print("# x0 -> x{}\n", i+1));
    M = Intersection(*M, *N);
  }

  Equality(*M, {0, k});

  // make sure cycle is a proper k-cycle
  for (auto i = 1U; i < k; i++) {
    if ((k % i) == 0) {
      dbg(OutputType::General, std::print("# x0 != x{}\n", i));
      M = Inequality(*M, {0, i});
    }
  }

  auto result = !M->Empty();
  dbg(OutputType::Quiet, std::print("Cycle({}, {}): {}\n", rule, k, result));
  return result;
}

// Construct a Büchi automaton to check if an ECA has a k-predecessor.
//
// A one-way infinite elementary cellular automaton has a k-predecessor if and
// only if there exist k distinct configurations x_1, ..., x_k evolving to a
// single configuration y after one application of the global map.
void Predecessor(int_type rule, int_type k, const std::vector<std::string>& p) {
  dbg(OutputType::General, std::print("# x0 -> y\n"));
  auto M = GlobalMap(rule, k+1, {0, k});

  std::unique_ptr<BuchiAutomaton> N;
  for (auto i = 1U; i < k; i++) {
    dbg(OutputType::General, std::print("# x{} -> y\n", i));
    N = GlobalMap(rule, k+1, {i, k});

    dbg(OutputType::General, std::print("# x0, ..., x{} -> y\n", i));
    M = Intersection(*M, *N);
  }

  for (auto i = 0UL; i < k-1; i++) {
    for (auto j = i+1; j < k; j++) {
      dbg(OutputType::General, std::print("# [x -> y] && [x{} != x{}]\n", i, j));
      M = Inequality(*M, {i, j});
    }
  }

  std::print("{:d}  ", !M->Empty());

  // ensure target isn't forbidden
  for (auto i = 0UL; i < p.size(); i++) {
    dbg(OutputType::General, std::print("# y != {}\n", p[i]));
    N = Pattern(k+1, k, p[i]);

    dbg(OutputType::General, std::print("# [x -> y] && [y != {}]\n", p[i]));
    M = Intersection(*M, *N);

    std::print("{:d}  ", !M->Empty());
  }
}

// A rule r is nilpotent after k applications if and only if there exists a configuration x_0 such that for all configurations u, u evolves to the same configuration y after k applications of the global map, and y is a fixed point
bool Nilpotent(int_type rule, int_type k) {
  BOOST_ASSERT_MSG(rule <= 255, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k >= 1, "Parameter k must be at least 1");

  // RuleTable(rule);

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, std::print("# Nilpotent({}, {})\n", rule, k));
  dbg(OutputType::General, std::print("# x0 -> x1\n"));
  auto M = GlobalMap(rule, k+2, {0, 1}, opts);

  for (auto i = 1U; i <= k; i++) {
    dbg(OutputType::General, std::print("# x{} -> x{}\n", i, i+1));
    auto N = GlobalMap(rule, k+2, {i, i+1}, opts);

    dbg(OutputType::General, std::print("# x0 -> x{}\n", i+1));
    M = Intersection(*M, *N);
  }

  // make sure that x_k is a fixed point
  dbg(OutputType::General, std::print("# [x0 -> x{}] && [x{} == x{}]\n", k+1, k, k+1));
  Equality(*M, {k, k+1});

  for (auto i = 0U; i <= k; i++) {
    M->ProjectLabel();
  }

  // M->final_states.reset(0);

  dbg(OutputType::General, M->Print());

  if (M->num_vertices == 0) {
    dbg(OutputType::Quiet, std::print("Nilpotent({}, {}): false\n", rule, k));
    return false;
  }

  auto map = M->Map();
  auto R = M->Determinize(*map);

  dbg(OutputType::General, R->Print());

  // R->Clean();
  // R->Minimize();
  // R->Clean();

  auto result = R->Universal();
  dbg(OutputType::Quiet, std::print("Nilpotent({}, {}): {}\n", rule, k, result));

  return result;
}

// Construct a Büchi automaton to check if an ECA has in-degree k
//
// A one-way infinite elementary cellular automaton has in-degree k if and
// only if there exist k distinct configurations x_1, ..., x_k evolving to a
// single configuration y after one application of the global map, and any
// other configuration u does not evolve to y

// \forall u, \exists x_1, ..., x_k, y:
// [(x_1 -> y) && ... && (x_k -> y) && (x_1 != x_2) && ... && (x_{k-1} != x_k)
//   && ((u -> y) => (u = x_1) || ... || (u = x_k))]
// Equivalently:
// \not\exists u, \exists x_1, ..., x_k, y:
// [(x_1 -/> y) || ... || (x_k -/> y) || (x_1 == x_2) || ... || (x_{k-1} == x_k)
//   || ((u -> y) && (u != x_1) && ... && (u != x_k))]
bool InDegree(int_type rule, int_type k) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k > 0, "Parameter k must be at least 1");

  // Track  Symbol
  // 0      y
  // 1      x_1
  // ...    ...
  // k      x_k
  // k+1    u

  dbg(OutputType::General, std::print("# InDegree({}, {})\n", rule, k));

  GlobalMapOpts opts;
  opts.negated = true;

  // The u term is (u -> y) && (u != x_1) && ... && (u != x_k), which is the
  // negation of the implication ((u -> y) => (u == x_1) || ... || (u == x_k)).
  // It uses the positive global map, not the negated one.
  dbg(OutputType::General, std::print("# u -> y\n"));
  auto M = GlobalMap(rule, k+2, {k+1, 0});

  for (auto i = 1UL; i <= k; i++) {
    dbg(OutputType::General, std::print("# [u -> y] && (x{} != u)\n", i));
    M = Inequality(*M, {i, k+1});
  }

  for (auto i = 1U; i <= k; i++) {
    dbg(OutputType::General, std::print("# x{} -/> y\n", i));
    auto N = GlobalMap(rule, k+2, {i, 0}, opts);

    dbg(OutputType::General, std::print("# [u -> y] || (x{} -/> y)\n", i));
    M = DisjointUnion(*M, *N);
  }

  for (auto i = 1U; i < k; i++) {
    for (auto j = i+1; j <= k; j++) {
      dbg(OutputType::General, std::print("# x{} == x{}\n", i, j));
      auto N = Equality(k+2, {i, j});
      dbg(OutputType::General, std::print("# [u -> y] || [x -/> y] || (x{} == x{})\n", i, j));
      M = DisjointUnion(*M, *N);
    }
  }

  dbg(OutputType::General, std::print("# [u -> y] || [x -/> y] || [x_i == x_j]\n"));

  // Remove track k+1 (u)
  M->ProjectLabel();

  auto map = M->Map();
  auto R = M->Determinize(*map);

  dbg(OutputType::General, R->Print());

  // R->Minimize();
  // R->Clean();
  // R->Minimize();

  bool result = !R->Universal();
  dbg(OutputType::General, std::print("InDegree({}, {}): {}\n", rule, k, result));
  return result;
  // boost::write_graphviz(std::cout, R->graph, boost::default_writer(), make_label_writer(boost::get(boost::edge_name, R->graph)));
}

// The global map is injective for a given rule r iff:
//   \forall    x,y,z: [(x -> z) && (y -> z) => (x == y)]
//   \notexists x,y,z: [(x -> z) && (y -> z) && (x != y)]
bool Injective(int_type rule) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, std::print("# Injective({})\n", rule));
  dbg(OutputType::General, std::print("# x -> z\n"));
  auto M = GlobalMap(rule, 3, {0, 2}, opts);

  dbg(OutputType::General, std::print("# y -> z\n"));
  auto N = GlobalMap(rule, 3, {1, 2}, opts);

  dbg(OutputType::General, std::print("# (x -> z) && (y -> z)\n"));
  M = Intersection(*M, *N);

  dbg(OutputType::General, std::print("# (x -> z) && (y -> z) && (x != y)\n"));
  M = Inequality(*M, {0, 1});

  // M->Minimize();

  auto result = M->Empty();
  dbg(OutputType::Quiet, std::print("Injective({}): {}\n", rule, result));
  return result;
}

// The global map for a given rule r is surjective iff:
//   \forall y, \exists x : x -> y
bool Surjective(int_type rule) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, std::print("# Surjective({})\n", rule));
  dbg(OutputType::General, std::print("# x -> y\n"));
  auto M = GlobalMap(rule, 2, {1, 0}, opts);

  dbg(OutputType::General, std::print("# _ -> y\n"));
  M->ProjectLabel();

  auto map = M->Map();
  auto R = M->Determinize(*map);

  dbg(OutputType::General, R->Print());

  // R->Clean();
  // R->Minimize();

  auto result = R->Universal();
  dbg(OutputType::Quiet, std::print("Surjective({}): {}\n", rule, result));
  return result;
}

// Generate the canonical one-step verifying automaton
void Canonical(int_type r) {
  BuchiAutomaton B(DE_BRUIJN_ALPHABET, DE_BRUIJN_SIZE);
  TransitionMap de_bruijn;

  RuleTable(r);

  std::print("# x -> y\n");
  auto M = GlobalMap(r, 2, {0, 1});

  // 8 transitions in total
  for (auto i = 0UL; i < DE_BRUIJN_TRANSITIONS; i++) {
    de_bruijn.insert({{source(i), map_bit(r, i)}, target(i, 1)});
  }

  B.Init(de_bruijn);
  B.initial_states.set(0);
  B.initial_states.set(1);
  B.final_states.set();

  std::print("# DE BRUIJN\n");
  B.Print();

  exit(EXIT_SUCCESS);
}

// Apply the global map of a given rule for a given pattern n of width k.
int_type Map(int_type rule, int_type n, int_type k) {
  constexpr static auto mask = 0x7U;

  if (k == 0) {
    return n;
  }

  auto width = 2*k-1;
  auto step = 0;
  auto count = 0UL;

  dbg(OutputType::Debug, std::print("# Map({}, {:0{}b}, {})\n", rule, n, 2*k+1, k));
  if (verbose > std::to_underlying(OutputType::Debug) && k > 1) {
    std::print("# p({:0{}b})\n", n, 2*k+1);
  }

  while (count < width) {
    dbg(OutputType::Debug, std::print("p({:0{}b})  =  {}\n", n & mask, 3, map_bit(rule, n & mask)););

    step = step | (map_bit(rule, n & mask) << count);

    n = n >> 1;
    count++;
  }

  dbg(OutputType::Debug, std::print("\n"));

  return Map(rule, step, k-1);
}

// Apply the transducer for the given rule to the language accepted by M.
std::unique_ptr<BuchiAutomaton> Cover(int_type rule, const BuchiAutomaton& M) {
  std::queue<int_pair> queue;

  // Q_B is the cartesian product of Q_M and {00,01,10,11}
  auto size = 4*M.num_vertices;

  boost::dynamic_bitset<> visited(size);

  auto B = std::make_unique<BuchiAutomaton>(2, size);
  TransitionMap de_bruijn;

  auto index = boost::get(boost::vertex_index, M.graph);
  auto label = boost::get(boost::edge_name, M.graph);

  for (auto i : dynamic_bitset_iterator(M.initial_states)) {
    auto u = boost::vertex(i, M.graph);

    for (auto [e_itr, e_end] = boost::out_edges(u, M.graph);
         e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, M.graph);
      auto q_M = index[v];

      auto a = label[*e_itr];

      queue.push({q_M, a});

      visited.set(compose_3(q_M, 0, a));
      B->initial_states.set(compose_3(q_M, 0, a));
    }
  }

  boost::dynamic_bitset<> S = ~M.final_states;
  B->final_states.set();

  for (auto i : dynamic_bitset_iterator(S)) {
    B->final_states.reset(compose_3(i, 0, 0));
    B->final_states.reset(compose_3(i, 0, 1));
    B->final_states.reset(compose_3(i, 1, 0));
    B->final_states.reset(compose_3(i, 1, 1));
  }

  while (!queue.empty()) {
    auto [q_M, q_T] = queue.front();

    queue.pop();

    dbg(OutputType::Debug, std::print("({}, {})\n", q_M, q_T));

    auto u = boost::vertex(q_M, M.graph);
    auto q = compose_3(q_M, 0, q_T);

    for (auto [e_itr, e_end] = boost::out_edges(u, M.graph); e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, M.graph);
      q_M = index[v];

      // M: p --> (a)   p'
      // T: q --> (a/b) q'
      //
      // M': (p, q) --> (b) (p', q')
      // The input is a, and the output is b = rho(q : a).

      auto a = label[*e_itr];
      auto b = Map(rule, compose_3(0, q_T, a), 1);
      auto k = compose_3(q_M, q_T & 0x1, a);

      // The transition leaves q_T and goes to (q_T & 0x1) : a
      de_bruijn.insert({{q, b}, k});

      dbg(OutputType::Debug, std::print("  -->  ({}, {})\n\n", q_M, compose_3(0, q_T & 0x1, a)));

      if (!visited[k]) {
        visited.set(k);
        queue.push({q_M, compose_3(0, q_T & 0x1, a)});
      }
    }
  }

  B->Init(de_bruijn);

  dbg(OutputType::Debug, B->Print());

  auto map = B->Map();
  auto N = B->RabinScott(*map);

  N->Minimize();

  return N;
}

// Creates a one-state Büchi automaton that accepts all configurations with an alphabet of size k.
std::unique_ptr<BuchiAutomaton> Top(int_type k) {
  TransitionMap de_bruijn;

  for (auto i = 0UL; i < k; i++) {
    de_bruijn.insert({{0, i}, 0});
  }

  auto B = std::make_unique<BuchiAutomaton>(k, 1);
  B->Init(de_bruijn);

  B->initial_states.set(0);
  B->final_states.set();

  dbg(OutputType::Debug, B->Print());

  return B;
}

// Generates the k-cover of a global map.
void Minimal(int_type rule, int_type k) {
  auto B = Top(2);

  for (auto i = 0U; i < k; i++) {
    B = Cover(rule, *B);
  }

  std::print("{},", B->num_vertices);
}

/**
  * Construct a Büchi automaton to verify a property of the specified one-way
  * infinite elementary cellular automaton.
 **/
void Run(int_type r, int_type k, const std::vector<std::string>& p) {
  dbg(OutputType::Debug, {
    std::print("rule = {}\n", r);
    std::print("k    = {}\n\n", k);
    std::print("p    = {{\n");
    for (auto& s : p) {
      std::print("  {}\n", s);
    }
    std::print("}}\n\n");

    RuleTable(r);
  });

  std::print("{:d}\n", Injective(r));
  std::print("{:d}\n", Surjective(r));
  std::print("{:d}\n", Cycle(r, k));
  // Predecessor(r, k, p);
  std::print("{:d}\n", InDegree(r, k));
  std::print("{:d}\n", Nilpotent(r, k));
  std::print("{:d}\n", Shift(r, k, ShiftType::Left));
  std::print("{:d}\n", Shift(r, k, ShiftType::Right));
  // Minimal(r, k);
  // Cover(r, nullptr);
  std::print("\n");
}

void Tabulate(int_type k) {
  std::vector<std::string> f0 = { "f0" };
  std::vector<std::string> f1 = { "f1" };
  std::vector<std::string> x;
  std::vector<std::string> y = { "0", "1" };
  std::vector<std::string> z =
    { "0", "1", "01", "001", "011", "0001", "0011", "0111" };

  // z.emplace_back("00001");
  // z.emplace_back("00011");
  // z.emplace_back("00101");
  // z.emplace_back("01001");
  // z.emplace_back("00111");
  // z.emplace_back("01011");
  // z.emplace_back("01111");

  // for (auto i = 0; i < 256; i++) {
  //   std::print("{:3}  ", i);
  //   std::print("{} {} ", binary_digits(i), decimal_digits(i));
  //   std::print("\n");
  // }

  for (auto i = 0; i < 256; i++) {
    std::print("{:3}", i);
    std::print("  {:d}", Injective(i));
    std::print("  {:d}", Surjective(i));
    // std::print("  {:d}", FixedPoint(i, y));
    // std::print("  {:d}", FixedPoint(i, z));
    // std::print("  {:d}", FixedPoint(i, f0));
    // std::print("  {:d}", FixedPoint(i, f1));
    std::print("  {:d}", Cycle(i, 1));
    std::print("  {:d}", Cycle(i, 2));
    std::print("  {:d}", Cycle(i, 3));
    std::print("  {:d}", Cycle(i, 4));
    std::print("  {:d}", Cycle(i, 5));
    std::print("  {:d}", Cycle(i, 6));
    // std::print("  {:d}", Cycle(i, 7));
    // std::print("  {:d}", Cycle(i, 8));
    // std::print("  {:d}", Cycle(i, 9));
    std::print("  {:d}", Nilpotent(i, 1));
    std::print("  {:d}", Nilpotent(i, 2));
    std::print("  {:d}", Nilpotent(i, 3));
    std::print("  {:d}", Nilpotent(i, 4));
    // std::print("  {:d}", Nilpotent(i, 5));
    std::print("  {:d}", InDegree(i, 1));
    std::print("  {:d}", InDegree(i, 2));
    std::print("  {:d}", InDegree(i, 3));
    // Predecessor(i, 1);
    // Predecessor(i, 2);
    // Predecessor(i, 2, x);
    // Predecessor(i, 3, x);
    // Predecessor(i, 3, z);
    // Predecessor(i, 4, x);
    // Predecessor(i, 5, x);
    std::print("  {:d}", Shift(i, 1, ShiftType::Left));
    std::print("  {:d}", Shift(i, 2, ShiftType::Left));
    std::print("  {:d}", Shift(i, 3, ShiftType::Left));
    std::print("  {:d}", Shift(i, 1, ShiftType::Right));
    std::print("  {:d}", Shift(i, 2, ShiftType::Right));
    std::print("  {:d}", Shift(i, 3, ShiftType::Right));
    // std::print("  {:d}", Shift(i, 3, ShiftType::Right, y));
    // std::print("  {:d}", Shift(i, 4, ShiftType::Right, y));
    // std::print("  {:d}", Shift(i, 5, ShiftType::Right, y));
    // Minimal(i, 1);
    // Minimal(i, 2);
    // Minimal(i, 3);
    // Minimal(i, 4);
    std::print("\n");
  }
}

} // namespace omega
