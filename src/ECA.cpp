/**
  * <ECA.cpp>
  *   Functions to construct and test properties of Büchi automata.
 **/

#include "BuchiAutomaton.h"
#include "ECA.h"
#include "Util.h"

#include <cmath>
#include <cstdio>
#include <cassert>
#include <iostream>
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
void RuleTable(uint8_t rule) {
  printf("# RULE TABLE %u\n", rule);

  // 8 transitions in total
  for (auto i = 0U; i < DE_BRUIJN_TRANSITIONS; i++) {
    printf("%u%u%u = %u\n", MAP(i, 2), MAP(i, 1), MAP(i, 0), MAP(rule, i));
  }

  printf("\n");
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
std::unique_ptr<BuchiAutomaton> GlobalMap(uint8_t rule, uint32_t k, int_pair pair, GlobalMapOpts opts = GlobalMapOpts()) {
  // S = 2^k, N = 10
  auto num_states = ECA_SIZE;
  if (!opts.full) {
    num_states--;
  }
  auto M = std::make_unique<BuchiAutomaton>(1 << k, num_states);

  if (verbose > static_cast<int>(OutputType::General)) {
    printf("# GLOBAL MAP\n");
    if (opts.negated) {
      printf("x%llu -/> x%llu", pair.first, pair.second);
    } else {
      printf("x%llu --> x%llu", pair.first, pair.second);
    }
    printf("  (N = %d, k = %d, S = 2^k = %llu)\n\n", num_states, k, M->num_alphabet);
  }

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
  M->initial_states.set(0);
  type[initial] = NodeType::Initial;

  uint8_t x1, x2, x3, y2, y3;

  auto width = binary_digits(M->num_alphabet);

  // 4 "initial" states of the form 0:x2:y2
  for (auto i = 0U; i < M->num_alphabet; i++) {
    // decompose i into two bits x2:y2
    x2 = MAP(i, pair.first);
    y2 = MAP(i, pair.second);

    // The numbering of states corresponds to the binary representation of its
    // component bits, accounting for the initial state.
    auto v = boost::vertex(COMPOSE_3(0, x2, y2) + 1, M->graph);

    dbg(OutputType::Debug, printf("0  "));
    dbg(OutputType::Debug, print_binary(i, width));
    dbg(OutputType::Debug, printf("  -->  [0%d %d]\n", x2, y2));

    auto [e, b] = boost::add_edge(initial, v, M->graph);
    label[e] = i;
  }

  dbg(OutputType::Debug, printf("\n"));

  auto [itr, end] = boost::vertices(M->graph);
  ++itr;

  // 8 possible combinations of x1:x2:y2 in total
  for (auto i = 1U; i < ECA_SIZE-1; i++) {
    auto index = i-1;
    auto u = *itr;
    ++itr;

    // decompose i into three bits x1:x2:y2
    x1 = GET_BIT(index, 2);
    x2 = GET_BIT(index, 1);
    y2 = GET_BIT(index, 0);

    // 4 possible combinations of x3:y3
    for (auto j = 0U; j < M->num_alphabet; j++) {
      // decompose j into two bits x3:y3
      x3 = GET_BIT(j, pair.first);
      y3 = GET_BIT(j, pair.second);
      auto v = boost::vertex(COMPOSE_3(x2, x3, y3) + 1, M->graph);

      dbg(OutputType::Debug, printf("%u  [%d%d %d]  ", index, x1, x2, y2));
      dbg(OutputType::Debug, print_binary(j, width));

      // x1:x2:y2 -> x2:x3:y3
      // if and only if MAP(x1, x2, x3) == y2
      if (MAP(rule, COMPOSE_3(x1, x2, x3)) == y2) {
        dbg(OutputType::Debug, printf("  -->  %d  [%d%d %d]\n", COMPOSE_3(x2, x3, y3), x2, x3, y3));
        dbg(OutputType::Debug, printf("  %d%d%d --> %d\n", x1, x2, x3, y2));

        auto [e, added] = boost::add_edge(u, v, M->graph);
        label[e] = j;
      } else {
        dbg(OutputType::Debug, printf("  -->  X\n"));

        if (opts.full) {
          auto [e, added] = boost::add_edge(u, crash, M->graph);
          label[e] = j;
        }
      }
    }

    dbg(OutputType::Debug, printf("\n"));

    if (opts.negated) {
      type[u] = NodeType::None;
    } else {
      M->final_states.set(i);
      type[u] = NodeType::Final;
    }
  }

  M->num_vertices = boost::num_vertices(M->graph);
  M->num_edges = boost::num_edges(M->graph);

  dbg(OutputType::Debug, printf("\n"));
  dbg(OutputType::Debug, M->Print());

  if (!opts.full) {
    M->Clean();

    // Since initial states are also final.
    // TODO: remove if changed.
    M->final_states.reset(0);

    dbg(OutputType::General, M->Print());
  }

  return M;
}

// Construct a Büchi automaton to check if two tracks are equal.
// k is the number of tracks
// pair indicates which tracks should be checked for equality
std::unique_ptr<BuchiAutomaton> Equality(uint32_t k, int_pair pair) {
  auto M = std::make_unique<BuchiAutomaton>(1 << k, 1);

  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  auto u = boost::vertex(0, M->graph);
  type[u] = NodeType::Both;

  for (auto i = 0U; i < M->num_alphabet; i++) {
    if (MAP(i, pair.first) == MAP(i, pair.second)) {
      auto [e, added] = boost::add_edge(u, u, M->graph);
      BOOST_ASSERT(added);

      label[e] = i;
    }
  }

  // M->final_states.set(0);

  M->Resize();

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  return M;
}

// Modify a Büchi automaton to ensure that two tracks are equal.
void Equality(BuchiAutomaton& M, int_pair pair) {
  // auto N = std::make_unique<BuchiAutomaton>(M.num_alphabet, M.num_vertices);

  auto label = boost::get(boost::edge_name, M.graph);
  auto type = boost::get(boost::vertex_name, M.graph);

  // add crash state
  // auto crash = boost::add_vertex(M.graph);
  // type[crash] = NodeType::None;

  auto done = false;
  while (!done) {
    done = true;

    for (auto [e_itr, e_end] = boost::edges(M.graph); e_itr != e_end; ++e_itr) {
      auto u = boost::source(*e_itr, M.graph);
      auto v = boost::target(*e_itr, M.graph);
      auto symbol = boost::get(boost::edge_name, M.graph, *e_itr);

      // if (v == crash) {
      //   continue;
      // }

      if (MAP(symbol, pair.first) != MAP(symbol, pair.second)) {
        boost::remove_edge(*e_itr, M.graph);

        // auto [e, added] = boost::add_edge(u, crash, M.graph);
        // label[e] = symbol;

        done = false;
        break;
      }
    }
  }

  // auto [e_itr, e_end] = boost::edges(M.graph);
  // for (auto e_next = e_itr; e_itr != e_end; e_itr = e_next) {
  //   ++e_next;

  //   auto u = boost::source(*e_itr, M.graph);
  //   auto v = boost::target(*e_itr, M.graph);
  //   auto symbol = boost::get(boost::edge_name, M.graph, *e_itr);

  //   if (v == crash) {
  //     continue;
  //   }

  //   if (MAP(symbol, pair.first) != MAP(symbol, pair.second)) {
  //     boost::remove_edge(*e_itr, M.graph);

  //     auto [e, b] = boost::add_edge(u, crash, M.graph);
  //     label[e] = symbol;
  //   }
  // }

  // for (auto i = 0U; i < M.num_alphabet; i++) {
  //   auto [e, added] = boost::add_edge(crash, crash, M.graph);
  //   label[e] = i;
  // }

  M.Resize();
  if (verbose > static_cast<int>(OutputType::General)) {
    M.Print();
  }

  M.Reachable();
  if (verbose > static_cast<int>(OutputType::General)) {
    M.Print();
  }
}

// Specialized product construction with the inequality automaton.
// The product P is two copies of M, one where the tracks have not differed,
// and one where they have.
std::unique_ptr<BuchiAutomaton> Inequality(BuchiAutomaton& M, int_pair pair) {
  auto P = std::make_unique<BuchiAutomaton>(M.num_alphabet);
  P->Reserve(2 * M.num_vertices);

  dbg(OutputType::Debug, printf("# (N = %llu, S = %llu)\n\n", M.num_vertices, M.num_alphabet));
  dbg(OutputType::Debug, printf("# Initial States\n"));

  auto vertex_width = decimal_digits(M.num_vertices);
  char fmt_vrtx[MAXLINE];
  snprintf(fmt_vrtx, MAXLINE, "(%%%du, %%d)\n", vertex_width);

  std::unordered_map<int_pair, graph_t::vertex_descriptor> map;
  std::queue<std::tuple<int_type, int_type, graph_t::vertex_descriptor>> queue;

  auto label = boost::get(boost::edge_name, P->graph);
  auto type = boost::get(boost::vertex_name, P->graph);

  auto label_M = boost::get(boost::edge_name, M.graph);
  auto index_M = boost::get(boost::vertex_index, M.graph);

  auto EQ = false;
  auto NE = true;

  ITERATE_BITSET(i, M.initial_states) {
    // Add a vertex to the graph representing the state (i, iEQ).
    auto u = boost::add_vertex(P->graph);

    // Store (i, EQ) -> u in the map.
    map.insert({{i, EQ}, u});

    // Fill in the graph as a BFS from the initial vertices.
    queue.push({i, EQ, u});

    // The state (i, EQ) is an initial state in P.
    type[u] = NodeType::Initial;

    dbg(OutputType::Debug, printf(fmt_vrtx, i, EQ));
  }

  dbg(OutputType::Debug, printf("\n"));

  char fmt_edge[MAXLINE];
  snprintf(fmt_edge, MAXLINE, "  -->  (%%%du, %%d)", vertex_width);

  auto edge_width = binary_digits(P->num_alphabet);

  while (!queue.empty()) {
    auto [source, source_component, u] = queue.front();
    queue.pop();

    dbg(OutputType::Debug, printf(fmt_vrtx, source, source_component));

    for (auto [e_itr, e_end] = boost::out_edges(boost::vertex(source, M.graph), M.graph); e_itr != e_end; ++e_itr) {
      auto symbol = label_M[*e_itr];

      auto target = index_M[boost::target(*e_itr, M.graph)];
      auto component = source_component;

      // If the tracks are not equal, go into the second component.
      if (GET_BIT(symbol, pair.first) != GET_BIT(symbol, pair.second)) {
        component = NE;
      }

      if (verbose > static_cast<int>(OutputType::Debug)) {
        printf("  ");
        print_binary(symbol, edge_width);
        printf(fmt_edge, target, component);
      }

      auto itr = map.find({target, component});

      // Add a new vertex to the product machine.
      if (itr == map.end()) {
        auto v = boost::add_vertex(P->graph);

        if (component == NE && M.final_states.test(target)) {
          type[v] = NodeType::Final;
        } else {
          type[v] = NodeType::None;
        }

        dbg(OutputType::Debug, printf("  *"));
        itr = map.insert({{target, component}, v}).first;

        queue.push({target, component, v});
      }

      // Add a new edge from u to the target vertex.
      auto [e, added] = boost::add_edge(u, itr->second, P->graph);
      label[e] = symbol;

      dbg(OutputType::Debug, printf("\n"));
    }

    dbg(OutputType::Debug, printf("\n"));
  }

  P->Resize();

  // Remove all initial states from the set of final states.
  // P->final_states &= ~P->initial_states;

  dbg(OutputType::General, P->Print());

  return P;
}

// Construct a Büchi automaton to check if two tracks are not equal.
std::unique_ptr<BuchiAutomaton> Inequality(uint32_t k, int_pair pair) {
  auto M = std::make_unique<BuchiAutomaton>(1 << k, 2);

  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  auto u = boost::vertex(0, M->graph);
  auto v = boost::vertex(1, M->graph);

  type[u] = NodeType::Initial;
  type[v] = NodeType::Final;

  for (auto i = 0U; i < M->num_alphabet; i++) {
    if (GET_BIT(i, pair.first) == GET_BIT(i, pair.second)) {
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

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  return M;
}

// Construct a Büchi automaton to check if track is a forbidden configuration.
std::unique_ptr<BuchiAutomaton> Pattern(uint32_t k, uint32_t n, const std::string& s) {
  auto start = 0U;

  graph_t::vertex_descriptor u, v, w;
  graph_t::edge_descriptor e;
  bool added;

  auto M = std::make_unique<BuchiAutomaton>(1 << k, s.length() + 2);
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
      uint32_t symbol = s[i % s.length()] - '0';
      if (MAP(j, n) == symbol) {
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

  uint32_t symbol;
  for (auto j = 0UL; j < M->num_alphabet; j++) {
    if (s[start + 1] != '*') {
      symbol = s[start] - '0';
      if (MAP(j, n) == symbol) {
        auto [e, b] = boost::add_edge(v, w, M->graph);
      } else if (s[s.length() - 1] != '*') {
        auto [e, b] = boost::add_edge(v, u, M->graph);
      }

      label[e] = j;
    } else {
      symbol = s[0] - '0';
      if (MAP(j, n) == symbol) {
        auto [e, b] = boost::add_edge(v, v, M->graph);
        label[e] = j;
      }

      if (MAP(j, n) == (uint32_t) (s[2] - '0')) {
        auto [e, b] = boost::add_edge(v, w, M->graph);
        label[e] = j;
      } else if (MAP(j, n) != symbol) {
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

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  M->Clean();

  if (s[0] != '2') {
    M->final_states.reset(0);
  }

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  return M;
}

// Construct a Büchi automaton that checks if a track is a finite configuration (but has at least one).
std::unique_ptr<BuchiAutomaton> Finite(uint32_t k, uint32_t n, const std::string& s) {
  auto num_states = 2;
  auto M = std::make_unique<BuchiAutomaton>(1 << k, num_states);
  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  auto u = boost::vertex(0, M->graph);
  type[u] = NodeType::Initial;

  auto v = boost::vertex(1, M->graph);
  type[v] = NodeType::Final;

  // finite occurrence of s[1]
  for (auto j = 0UL; j < M->num_alphabet; j++) {
    if (MAP(j, n) != (uint32_t) (s[1] - '0')) {
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

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  auto N = std::make_unique<BuchiAutomaton>(1 << k, 2);
  label = boost::get(boost::edge_name, N->graph);
  type = boost::get(boost::vertex_name, N->graph);

  u = boost::vertex(0, N->graph);
  type[u] = NodeType::Initial;

  v = boost::vertex(1, N->graph);
  type[v] = NodeType::Final;

  // at least one occurrence of s[1]
  for (auto j = 0UL; j < N->num_alphabet; j++) {
    if (MAP(j, n) != (uint32_t) (s[1] - '0')) {
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

// Construct a Büchi automaton that checks if track1 is a right shift of track2.
std::unique_ptr<BuchiAutomaton> RightShift(uint32_t k, int_pair pair, bool partial = true) {
  auto num_states = 4;
  if (partial) {
    num_states--;
  }

  auto M = std::make_unique<BuchiAutomaton>(1 << k, num_states);
  auto label = boost::get(boost::edge_name, M->graph);
  auto type = boost::get(boost::vertex_name, M->graph);

  // initial state
  auto initial = boost::vertex(0, M->graph);
  type[initial] = NodeType::Initial;

  // crash state
  graph_t::vertex_descriptor x;
  if (!partial) {
    x = boost::vertex(3, M->graph);
    type[x] = NodeType::None;
  }

  for (auto i = 0UL; i < M->num_alphabet; i++) {
    auto n = MAP(i, pair.first);
    auto u = boost::vertex(n+1, M->graph);

    // edge from initial state
    auto [e, added] = boost::add_edge(initial, u, M->graph);
    label[e] = i;

    // edge on crash state
    if (!partial) {
      auto [e, added] = boost::add_edge(x, x, M->graph);
      label[e] = i;
    }
  }

  for (auto i = 1U; i <= 2; i++) {
    auto u = boost::vertex(i, M->graph);
    type[u] = NodeType::Final;

    for (auto j = 0U; j < M->num_alphabet; j++) {
      auto n = MAP(j, pair.first);
      auto v = boost::vertex(n+1, M->graph);

      // previous symbol of x matches the current symbol of y
      if (MAP(j, pair.second) == i) {
        auto [e, added] = boost::add_edge(u, v, M->graph);
        label[e] = j;
      } else if (!partial) {
        auto [e, added] = boost::add_edge(u, x, M->graph);
        label[e] = j;
      }
    }
  }

  M->Resize();
  // M->initial_states.set(0);

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  return M;
}

bool RightShift(uint32_t rule, uint32_t k, const std::vector<std::string>& p) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k >= 1, "Parameter k must be at least 1");

  GlobalMapOpts opts;
  opts.full = false;
  auto M = GlobalMap(rule, k+1, {0, 1}, opts);

  for (auto i = 1U; i < k; i++) {
    auto N = GlobalMap(rule, k+1, {i, i+1}, opts);
    M = Intersection(*M, *N);
  }

  auto N = RightShift(k+1, {0, k});
  M = Intersection(*M, *N);

  // printf("%d  ", !M->Empty());

  // for (const auto& s : p) {
  //   dbg(OutputType::General, printf("# x != %s\n", s.c_str()));
  //   N = Pattern(k+1, 0, s);

  //   dbg(OutputType::General, printf("# [x -> y] && x != %s\n", s.c_str()));
  //   M = Intersection(*M, *N);

  //   printf("%d  ", !M->Empty());
  // }

  auto result = !M->Empty();
  if (verbose > static_cast<int>(OutputType::General)) {
    std::cout << "RightShift(" << rule << ", " << k << "): " << std::boolalpha << result << std::endl;
  }
  return result;
}

// Construct a Büchi automaton to check if an ECA has a fixed-point.
// A one-way infinite elementary cellular automaton has a fixed-point if there exists a configuration that evolves to itself after one application of the global map.
bool FixedPoint(uint32_t rule) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, printf("# Fixed Point %u\n", rule));

  // Construct x_0 -> x_1 and x_0 == x_1
  dbg(OutputType::General, printf("# x0 -> x1\n"));
  auto M = GlobalMap(rule, 2, {0, 1}, opts);

  dbg(OutputType::General, printf("# x0 == x1\n"));
  Equality(*M, {0, 1});

  auto result = !M->Empty();
  if (verbose > static_cast<int>(OutputType::Quiet)) {
    std::cout << "FixedPoint(" << rule << "): " << std::boolalpha << result << std::endl;
  }
  return result;

  // std::unique_ptr<BuchiAutomaton> N;
  // ensure fixed point isn't forbidden
  // for (auto i = 0UL; i < p.size(); i++) {
  //   // finite support
  //   if (p[i][0] == 'f') {
  //     if ((p[i][1] == '1') && MAP(rule, 0)) {
  //       printf("X  ");
  //       continue;
  //     }

  //     else if ((p[i][1] == '0') && !MAP(rule, 7)) {
  //       printf("X  ");
  //       continue;
  //     }

  //     dbg(OutputType::General, printf("# x0 == %s\n", p[i].c_str()));
  //     N = Finite(2, 0, p[i]);
  //     dbg(OutputType::General, printf("# [x -> x] && x0 == %s\n", p[i].c_str()));
  //   } else {
  //     dbg(OutputType::General, printf("# x0 != %s\n", p[i].c_str()));
  //     N = Pattern(2, 0, p[i]);
  //     dbg(OutputType::General, printf("# [x -> x] && x0 != %s\n", p[i].c_str()));
  //   }

  //   M = Intersection(*M, *N);

  //   // printf("%d  ", !M->Empty());
  //   return !M->Empty();
  // }
}

// Construct a Büchi automaton to check if an ECA has a k-cycle.
//
// A one-way infinite elementary cellular automaton has a k-cycle if there exists a configuration that evolves to itself after k applications of the global map.
//
// Furthermore, the cycle must be proper, i.e., not a d-cycle where d is a divisor of k.
bool Cycle(uint32_t rule, uint32_t k) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k >= 1, "Parameter k must be at least 1");

  GlobalMapOpts opts;
  opts.full = false;

  // construct x_1 -> x_2
  // for i = {2, ..., k-1}, construct x_i -> x_{i+1} and x_1 -> x_{i+1}
  // finally, construct x_k -> x_1 and x_1 -> x_1 from x_1 -> x_k

  dbg(OutputType::General, printf("# Cycle(%d, %d)\n", rule, k));
  dbg(OutputType::General, printf("# x0 -> x1\n"));
  auto M = GlobalMap(rule, k+1, {0, 1}, opts);

  for (auto i = 1U; i < k; i++) {
    dbg(OutputType::General, printf("# x%u -> x%u\n", i, i+1));
    auto N = GlobalMap(rule, k+1, {i, i+1}, opts);

    dbg(OutputType::General, printf("# x0 -> x%u\n", i+1));
    M = Intersection(*M, *N);
  }

  dbg(OutputType::General, printf("# x0 == x%u\n", k));
  Equality(*M, {0, k});

  // make sure cycle is a proper k-cycle
  for (auto i = 1U; i < k; i++) {
    if ((k % i) == 0) {
      dbg(OutputType::General, printf("# x0 != x%u\n", i));
      M = Inequality(*M, {0, i});
    }
  }

  auto result = !M->Empty();
  if (verbose > static_cast<int>(OutputType::Quiet)) {
    std::cout << "Cycle(" << rule << ", " << k << "): " << std::boolalpha << result << std::endl;
  }
  return result;
}

// Construct a Büchi automaton to check if an ECA has a k-predecessor.
//
// A one-way infinite elementary cellular automaton has a k-predecessor if and
// only if there exist k distinct configurations x_1, ..., x_k evolving to a
// single configuration y after one application of the global map.
void Predecessor(uint8_t rule, uint32_t k, const std::vector<std::string>& p) {
  dbg(OutputType::General, printf("# x0 -> y\n"));
  auto M = GlobalMap(rule, k+1, {0, k});

  std::unique_ptr<BuchiAutomaton> N;
  for (auto i = 1U; i < k; i++) {
    dbg(OutputType::General, printf("# x%u -> y\n", i));
    N = GlobalMap(rule, k+1, {i, k});

    dbg(OutputType::General, printf("# x0, ..., x%u -> y\n", i));
    M = Intersection(*M, *N);
  }

  for (auto i = 0UL; i < k-1; i++) {
    for (auto j = i+1; j < k; j++) {
      dbg(OutputType::General, printf("# [x -> y] && [x%zd != x%zd]\n", i, j));
      M = Inequality(*M, {i, j});
    }
  }

  printf("%d  ", !M->Empty());

  // ensure target isn't forbidden
  for (auto i = 0UL; i < p.size(); i++) {
    dbg(OutputType::General, printf("# y != %s\n", p[i].c_str()));
    N = Pattern(k+1, k, p[i]);

    dbg(OutputType::General, printf("# [x -> y] && [y != %s]\n", p[i].c_str()));
    M = Intersection(*M, *N);

    printf("%d  ", !M->Empty());
  }
}

bool Nilpotent(uint32_t rule, uint32_t k) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k > 0, "Parameter k must be at least 1");

  GlobalMapOpts opts;
  opts.full = false;

  // construct x_1 -> x_2
  // for i = {2, ..., k+1}, construct x_i -> x_{i+1} and x_1 -> x_{i+1}
  // at the end, the machine is x_1 -> x_{k+1}
  // dbg(OutputType::General, printf("# x%u -> x%u\n", k, k+1));
  // auto M = GlobalMap(rule, k+2, {k, k+1}, opts);

  // // make sure that x_k is a fixed point
  // dbg(OutputType::General, printf("# [x%u -> x%u] && [x%u == x%u]\n", k, k+1, k, k+1));
  // Equality(*M, {k, k+1});

  // for (auto i = 0U; i < k; i++) {
  //   dbg(OutputType::General, printf("# x%u -> x%u\n", i, i+1));
  //   auto N = GlobalMap(rule, k+2, {i, i+1}, opts);

  //   dbg(OutputType::General, printf("# x0 -> x%u\n", i+1));
  //   M = Intersection(*M, *N);
  // }

  dbg(OutputType::General, printf("# x0 -> x1\n"));
  auto M = GlobalMap(rule, k+1, {0, 1}, opts);

  for (auto i = 1U; i < k; i++) {
    dbg(OutputType::General, printf("# x%u -> x%u\n", i, i+1));
    auto N = GlobalMap(rule, k+1, {i, i+1}, opts);

    dbg(OutputType::General, printf("# x0 -> x%u\n", i+1));
    M = Intersection(*M, *N);
  }

  // make sure that x_k is a fixed point
  dbg(OutputType::General, printf("# [x0 -> x%u] && [x%u == x%u]\n", k, k-1, k));
  Equality(*M, {k-1, k});

  for (auto i = 0U; i < k; i++) {
    M->ProjectLabel();
  }

  M->final_states.reset(0);

  if (verbose > static_cast<int>(OutputType::General)) {
    M->Print();
  }

  auto map = M->Map();
  auto R = M->Determinize(*map);

  if (verbose > static_cast<int>(OutputType::General)) {
    R->Print();
  }

  // R->Clean();
  // R->Minimize();
  // R->Clean();

  auto result = R->Universal();
  if (verbose > static_cast<int>(OutputType::Quiet)) {
    std::cout << "Nilpotent(" << rule << ", " << k << "): " << std::boolalpha << result << std::endl;
  }

  return result;
}

// Construct a Büchi automaton to check if an ECA has in-degree k.
//
// A one-way infinite elementary cellular automaton has in-degree k if and
// only if there exist k distinct configurations x_1, ..., x_k evolving to a
// single configuration y after one application of the global map, and any
// other configuration u does not evolve to y.

// \exists x_1, ..., x_k, y, \forall u
// [(x_1 -> y) && ... && (x_k -> y) && (x_1 != x_2) && ... && (x_{k-1} != x_k)
//   && ((u -> y) => (u = x_1) || ... || (u = x_k))]

// \exists x_1, ..., x_k, y, \not\exists u
// [(x_1 -/> y) || ... || (x_k -/> y) || (x_1 == x_2) || ... || (x_{k-1} == x_k)
//   || ((u -/> y) && (u != x_1) && ... && (u != x_k))]
bool InDegree(uint32_t rule, uint32_t k) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");
  BOOST_ASSERT_MSG(k > 0, "Parameter k must be at least 1");

  // Track  Symbol
  // 0      y
  // 1      x_1
  // ...    ...
  // k      x_k
  // k+1    u

  GlobalMapOpts opts;
  // opts.full = false;
  opts.negated = true;
  dbg(OutputType::General, printf("# InDegree %d %d\n", rule, k));
  dbg(OutputType::General, printf("# u -/> y\n"));
  auto M = GlobalMap(rule, k+2, {k+1, 0}, opts);

  for (auto i = 1UL; i <= k; i++) {
    dbg(OutputType::General, printf("# [u -/> y] && (x%zd != u)\n", i));
    M = Inequality(*M, {i, k+1});
  }

  for (auto i = 1U; i <= k; i++) {
    dbg(OutputType::General, printf("# x%u -/> y\n", i));
    auto N = GlobalMap(rule, k+2, {i, 0}, opts);

    dbg(OutputType::General, printf("# [u -/> y] || (x%u -/> y)\n", i));
    M = DisjointUnion(*M, *N);
  }

  // for (auto i = 1UL; i <= k; i++) {
  //   dbg(OutputType::General, printf("# x%zd -/> y\n", i));
  //   auto N = GlobalMap(r, k+2, {i, 0}, opts);
  //   dbg(OutputType::General, printf("# [u -/> y] || (x%zd -/> y)\n", i));
  //   M = DisjointUnion(*M, *N);
  // }

  for (auto i = 1U; i < k; i++) {
    for (auto j = i+1; j <= k; j++) {
      dbg(OutputType::General, printf("# x%u == x%u\n", i, j));
      auto N = Equality(k+2, {i, j});
      dbg(OutputType::General, printf("# [u -/> y] || [x -/> y] || (x%u == x%u)\n", i, j));
      M = DisjointUnion(*M, *N);
    }
  }

  dbg(OutputType::General, printf("# [u -/> y] || [x -/> y] || [x_i == x_j]\n"));

  // Remove track k+1 (u).
  M->ProjectLabel();

  auto map = M->Map();
  auto R = M->Determinize(*map);

  dbg(OutputType::General, R->Print());

  // R->Minimize();
  // R->Clean();
  // R->Minimize();

  bool result = !R->Universal();
  if (verbose > static_cast<int>(OutputType::General)) {
    std::cout << "InDegree(" << rule << ", " << k << "): " << std::boolalpha << result << std::endl;
  }
  return result;
  // boost::write_graphviz(std::cout, R->graph, boost::default_writer(), make_label_writer(boost::get(boost::edge_name, R->graph)));
}

// The global map is injective for a given rule r iff:
//   \forall    x,y,z: [(x -> z) && (y -> z) => (x == y)]
//   \notexists x,y,z: [(x -> z) && (y -> z) && (x != y)]
bool Injective(uint32_t rule) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, printf("# x -> z\n"));
  auto M = GlobalMap(rule, 3, {0, 2}, opts);

  dbg(OutputType::General, printf("# y -> z\n"));
  auto N = GlobalMap(rule, 3, {1, 2}, opts);

  dbg(OutputType::General, printf("# (x -> z) && (y -> z)\n"));
  M = Intersection(*M, *N);

  dbg(OutputType::General, printf("# (x -> z) && (y -> z) && (x != y)\n"));
  M = Inequality(*M, {0, 1});

  // M->Minimize();

  auto result = M->Empty();
  if (verbose > static_cast<int>(OutputType::Quiet)) {
    std::cout << "Injective(" << rule << "): " << std::boolalpha << result << std::endl;
  }
  return result;
}

// The global map for a given rule r is surjective iff:
//   \forall y, \exists x : x -> y
bool Surjective(uint32_t rule) {
  BOOST_ASSERT_MSG(rule < 256, "Rule must be in the range [0, 255]");

  GlobalMapOpts opts;
  opts.full = false;

  dbg(OutputType::General, printf("# x -> y\n"));
  auto M = GlobalMap(rule, 2, {1, 0}, opts);

  dbg(OutputType::General, printf("# _ -> y\n"));
  M->ProjectLabel();

  auto map = M->Map();
  auto R = M->Determinize(*map);

  if (verbose > static_cast<int>(OutputType::General)) {
    R->Print();
  }

  // R->Clean();
  // R->Minimize();

  auto result = R->Universal();
  if (verbose > static_cast<int>(OutputType::Quiet)) {
    std::cout << "Surjective(" << rule << "): " << std::boolalpha << result << std::endl;
  }
  return result;
}

/**
  * Generate the canonical one-step verifying automaton.
 **/
void Canonical(uint8_t r) {
  BuchiAutomaton B(DE_BRUIJN_ALPHABET, DE_BRUIJN_SIZE);
  TransitionMap de_bruijn;

  RuleTable(r);

  printf("# x -> y\n");
  auto M = GlobalMap(r, 2, {0, 1});

  // 8 transitions in total
  for (auto i = 0UL; i < DE_BRUIJN_TRANSITIONS; i++) {
    de_bruijn.insert({{SOURCE(i), MAP(r, i)}, TARGET(i, 1)});
  }

  B.Init(de_bruijn);
  B.initial_states.set(0);
  B.initial_states.set(1);
  B.final_states.set();

  printf("# DE BRUIJN\n");
  B.Print();
  std::flush(std::cout);

  exit(EXIT_SUCCESS);
}

// Apply the global map of a given rule for a given pattern n of width k.
uint8_t Map(uint8_t rule, uint32_t n, uint32_t k) {
  const static auto mask = 0x7U;

  if (k == 0) {
    return n;
  }

  uint32_t width = 2*k-1;
  uint32_t step = 0;
  uint32_t count = 0;

  if (verbose > static_cast<int>(OutputType::Debug) && k > 1) {
    printf("p(");
    print_binary(n, 2*k+1);
    printf(")\n\n");
  }

  while (count < width) {
    if (verbose > static_cast<int>(OutputType::Debug)) {
      printf("p(");
      print_binary(n & mask, 3);
      printf(")  =  %u\n", MAP(rule, n & mask));
    }

    step = step | (MAP(rule, n & mask) << count);

    n = n >> 1;
    count++;
  }

  dbg(OutputType::Debug, std::cout << std::endl);

  return Map(rule, step, k-1);
}

// Apply the transducer for the given rule to the language accepted by M.
std::unique_ptr<BuchiAutomaton> Cover(uint8_t rule, const BuchiAutomaton& M) {
  std::queue<uint32_t> m_queue;
  std::queue<uint32_t> t_queue;

  // Q_B is the cartesian product of Q_M and {00,01,10,11}
  auto size = 4*M.num_vertices;

  boost::dynamic_bitset<> visited(size);

  auto B = std::make_unique<BuchiAutomaton>(2, size);
  TransitionMap de_bruijn;

  auto index = boost::get(boost::vertex_index, M.graph);
  auto label = boost::get(boost::edge_name, M.graph);

  ITERATE_BITSET(i, M.initial_states) {
    auto u = boost::vertex(i, M.graph);

    for (auto [e_itr, e_end] = boost::out_edges(u, M.graph);
         e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, M.graph);
      auto q_M = index[v];

      auto a = label[*e_itr];

      m_queue.push(q_M);
      t_queue.push(a);

      visited.set(COMPOSE_3(q_M, 0, a));
      B->initial_states.set(COMPOSE_3(q_M, 0, a));
    }
  }

  boost::dynamic_bitset<> S = ~M.final_states;
  B->final_states.set();

  ITERATE_BITSET(i, S) {
    B->final_states.reset(COMPOSE_3(i, 0, 0));
    B->final_states.reset(COMPOSE_3(i, 0, 1));
    B->final_states.reset(COMPOSE_3(i, 1, 0));
    B->final_states.reset(COMPOSE_3(i, 1, 1));
  }

  while (!m_queue.empty()) {
    auto q_M = m_queue.front();
    auto q_T = t_queue.front();

    m_queue.pop();
    t_queue.pop();

    dbg(OutputType::Debug, printf("(%d, %d)\n", q_M, q_T));

    auto u = boost::vertex(q_M, M.graph);
    auto q = COMPOSE_3(q_M, 0, q_T);

    for (auto [e_itr, e_end] = boost::out_edges(u, M.graph); e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, M.graph);
      q_M = index[v];

      // M: p --> (a)   p'
      // T: q --> (a/b) q'
      //
      // M': (p, q) --> (b) (p', q')
      // The input is a, and the output is b = rho(q : a).

      auto a = label[*e_itr];
      auto b = Map(rule, COMPOSE_3(0, q_T, a), 1);
      auto k = COMPOSE_3(q_M, q_T & 0x1, a);

      // The transition leaves q_T and goes to (q_T & 0x1) : a
      de_bruijn.insert({{q, b}, k});

      dbg(OutputType::Debug, printf("  -->  (%u, %u)\n\n", q_M, COMPOSE_3(0, q_T & 0x1, a)));

      if (!visited.test(k)) {
        visited.set(k);

        m_queue.push(q_M);
        t_queue.push(COMPOSE_3(0, q_T & 0x1, a));
      }
    }
  }

  // char fmt[MAXLINE];
  // snprintf(fmt, MAXLINE, "(%%%du, ", decimal_digits(M->num_vertices));

  B->Init(de_bruijn);

  if (verbose > static_cast<int>(OutputType::Debug)) {
    B->Print();
  }

  auto map = B->Map();
  auto N = B->RabinScott(*map);

  N->Minimize();

  return N;
}

// Creates a one-state Büchi automaton that accepts all configurations with an alphabet of size k.
std::unique_ptr<BuchiAutomaton> Top(uint8_t k) {
  TransitionMap de_bruijn;

  for (auto i = 0UL; i < k; i++) {
    de_bruijn.insert({{0, i}, 0});
  }

  auto B = std::make_unique<BuchiAutomaton>(k, 1);
  B->Init(de_bruijn);

  B->initial_states.set(0);
  B->final_states.set();

  if (verbose > static_cast<int>(OutputType::Debug)) {
    B->Print();
  }

  return B;
}

// Generates the k-cover of a global map.
void Minimal(uint8_t rule, uint32_t k) {
  auto B = Top(2);

  for (auto i = 0U; i < k; i++) {
    B = Cover(rule, *B);
  }

  printf("%llu,", B->num_vertices);
}

/**
  * Construct a Büchi automaton to verify a property of the specified one-way
  * infinite elementary cellular automaton.
 **/
void Run(uint8_t r, uint8_t k, const std::vector<std::string>& p) {
  if (verbose > static_cast<int>(OutputType::Debug)) {
    printf("rule = %d\n", r);
    printf("k    = %d\n\n", k);
    for (auto& s : p) {
      std::cout << s << "\n";
    }

    RuleTable(r);
  }

  printf("%d  ", Injective(r));
  printf("%d  ", Surjective(r));
  // FixedPoint(r, p);
  printf("%d  ", Cycle(r, k+1));
  // Predecessor(r, k, p);
  printf("%d  ", InDegree(r, k));
  // RightShift(r, k, p);
  printf("%d  ", Nilpotent(r, k));
  // Minimal(r, k);
  // Cover(r, nullptr);
  std::cout << std::endl;
}

void Tabulate(uint32_t k) {
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

  for (auto i = k - 1U; i < 256; i++) {
    printf("%3u  ", i);
    printf("%d  ", Injective(i));
    printf("%d  ", Surjective(i));
    printf("%d  ", Cycle(i, 2));
    printf("%d  ", Cycle(i, 3));
    // printf("%d  ", Cycle(i, 4));
    // printf("%d  ", Cycle(i, 5));
    // printf("%d  ", Cycle(i, 6));
    // printf("%d  ", Cycle(i, 7));
    // printf("%d  ", Cycle(i, 8));
    // printf("%d  ", Cycle(i, 9));
    printf("%d  ", Nilpotent(i, 1));
    // printf("%d  ", Nilpotent(i, 2));
    // printf("%d  ", Nilpotent(i, 3));
    // printf("%d  ", Nilpotent(i, 4));
    // printf("%d  ", Nilpotent(i, 5));
    printf("%d  ", InDegree(i, 1));
    printf("%d  ", InDegree(i, 2));
    // printf("%d  ", InDegree(i, 3));
    FixedPoint(i);
    // FixedPoint(i, y);
    // FixedPoint(i, z);
    // FixedPoint(i, f0);
    // FixedPoint(i, f1);
    // Predecessor(i, 1);
    // Predecessor(i, 2);
    // Predecessor(i, 2, x);
    // Predecessor(i, 3, x);
    // Predecessor(i, 3, z);
    // Predecessor(i, 4, x);
    // Predecessor(i, 5, x);
    RightShift(i, 2, x);
    // RightShift(i, 3, y);
    // RightShift(i, 4, y);
    // RightShift(i, 5, y);
    // Minimal(i, 1);
    // Minimal(i, 2);
    // Minimal(i, 3);
    // Minimal(i, 4);
    printf("\n");
  }
}

} // namespace omega
