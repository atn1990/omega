// RabinAutomaton.cpp
// *  Library to construct and manipulate Rabin automata.

#include "RabinAutomaton.h"
#include "Util.h"

#include <cstdio>
#include <iostream>

#include <boost/graph/strong_components.hpp>

template <>
struct std::formatter<omega::RabinPair> : std::formatter<std::string_view> {
    auto format(const omega::RabinPair& p, format_context& ctx) const {
        return std::format_to(ctx.out(), "L = {}\nR = {}\n", p.left, p.right);
    }
};

namespace omega {

std::ostream& operator<<(std::ostream& os, const RabinPair& p) {
  os << "L = " << p.left << "\nR = " << p.right << "\n\n";
  return os;
}

RabinAutomaton::RabinAutomaton(size_type alphabet, size_type vertices)
  : Automaton(alphabet, vertices) {}

// Build the underlying transition graph from a deterministic transition
// function represented by a map from vertex-symbol pairs to vertices.
void RabinAutomaton::Init(const RabinTransitionMap& map) {
  BOOST_ASSERT(map.size() == num_vertices * num_alphabet);

  const auto vertex_width = decimal_digits(num_vertices-1);
  const auto symbol_width = decimal_digits(num_alphabet-1);

  dbg(OutputType::Quiet, std::print("# TRANSITIONS\n"));

  auto label = boost::get(boost::edge_name, graph);

  // add edges in sorted order by source vertex and transition symbol
  for (auto i = 0UL; i < num_vertices; i++) {
    for (auto j = 0UL; j < num_alphabet; j++) {
      auto itr = map.find({i, j});

      // the transition system is deterministic
      BOOST_ASSERT(itr != map.end());
      dbg(OutputType::Quiet, std::print("({:{}}, {:{}})  -->  {:{}}\n", i, vertex_width, j, symbol_width, itr->second, vertex_width));

      auto [e, b] = boost::add_edge(i, itr->second, graph);
      label[e] = j;
    }
  }

  num_edges = boost::num_edges(graph);
}

void RabinAutomaton::Print() const {
  Automaton::Print();

  std::print("# PAIRS\n");
  for (const auto& pair : pairs) {
    std::print("{}\n", pair);
  }
}

/**
  * Tests an input string of the form uv* for acceptance, where u is an initial
  * transient segment and v is repeated infinitely often.
  *
  * The initial transient segment defines a path in the automaton starting from
  * the initial state. The repeating segment starts a walk in the graph until a
  * cycle is found. If the cycle intersects some Rabin condition, then the
  * input is accepted.
  *
  * TODO: Modify strings to be comma separated to accommodate larger alphabets.
 **/
void RabinAutomaton::TestUV(const std::string& U, const std::string& V) {
  // start with the initial state
  graph_t::vertex_descriptor v = boost::vertex(0, graph);
  graph_t::vertex_descriptor t;

  auto label = boost::get(boost::edge_name, graph);

  // iterate through the initial transient segment u
  for (auto c : U) {
    auto symbol = c - '0';

    // find edges out of current vertex and check if one matches input symbol
    for (auto [e_itr, e_end] = boost::out_edges(v, graph);
         e_itr != e_end; ++e_itr) {
      // if a matching edge is found, update current vertex to target of edge
      if (symbol == label[*e_itr]) {
        t = boost::target(*e_itr, graph);
        dbg(OutputType::Quiet, std::print("({}, {})  ->  {}\n", v, symbol, t));
        v = t;

        break;
      }
    }
  }

  size_t start = 0UL;

  // these arrays keep track of vertices and positions in V read in order
  // to determine if computation has entered a cycle
  std::vector<graph_t::vertex_descriptor> vertices;

  std::unordered_map<std::pair<graph_t::vertex_descriptor, size_t>, size_t> visited;

  // since V is repeated infinitely often, loop across V and reset i to the
  // beginning when the end is reached until a cycle is found

  bool has_cycle = false;
  for (auto i = 0UL; !has_cycle; i = (i + 1) % V.length()) {
    for (auto [e_itr, e_end] = boost::out_edges(v, graph);
         e_itr != e_end; ++e_itr) {
      uint32_t symbol = V[i] - '0';

      // when a matching edge is found, check if computation is in a cycle
      if (symbol != label[*e_itr]) {
        continue;
      }

      // keep track of every state seen so far as well as the position in V
      // when that state was seen

      // since computation has not reached a cycle, add the current vertex
      // and position in V to the list of previously visited states
      auto [itr, inserted] = visited.insert({{v, i}, vertices.size()});

      t = boost::target(*e_itr, graph);
      dbg(OutputType::Quiet, std::print("({}, {})  ->  {} ({})\n", v, symbol, t, i));
      v = t;

      if (inserted) {
        vertices.push_back(v);
        assert(
            vertices.size() <= (num_vertices * num_alphabet + 1) * V.length());
      } else {
        // If the pair (v, i) is not new, then the mapped value is the index
        // where it first appears in the vertex list.
        has_cycle = true;
        start = itr->second;
      }

      break;
    }
  }

  dbg(OutputType::Quiet, std::print("# Cycle found (transient = {})\n", vertices.size() - V.length()));

  // store cycle in a bitset for comparison against Rabin condition
  boost::dynamic_bitset<> cycle(num_vertices);
  auto index = boost::get(boost::vertex_index, graph);

  // the tail of the vertex list, beginning from start is part of a cycle
  for (auto i = start; i < vertices.size(); i++) {
    cycle.set(index[vertices[i]]);
    dbg(OutputType::Quiet, std::print("{}  ", vertices[i]));
  }

  dbg(OutputType::Quiet, std::print("\n{}\n", cycle));

  // once the cycle stored in a bitset, iterate through Rabin pairs and look
  // for L and R such that L is disjoint from C and R is non-disjoint from C
  for (const auto& pair : pairs) {
    if ((pair.left & cycle).none() && (pair.right & cycle).any()) {
      std::print("Input accepted\n");

      dbg(OutputType::Quiet, std::print("{}\n", pair));
      return;
    }
  }

  // if control reaches here, no accepting Rabin pair was found
  std::print("Input rejected\n");
}

/**
  * Removes useless or redundant Rabin pairs. A Rabin pair (L, R) may be
  * deleted if after removing all vertices of L, there is no non-trivial
  * strongly connected component intersecting R.
 **/
void RabinAutomaton::Clean() {
  dbg(OutputType::Quiet, std::print("# Clean\n\n"));

  const auto width = decimal_digits(num_vertices-1);
  auto pair = pairs.begin();
  while (pair != pairs.end()) {
    graph_t H(graph);
    auto index = boost::get(boost::vertex_index, H);

    // for each vertex in L, remove its corresponding edges
    for (auto i : dynamic_bitset_iterator(pair->left)) {
      boost::clear_vertex(boost::vertex(i, H), H);
    }

    dbg(OutputType::General, std::print("{}\n", *pair));

    // Compute the strongly connected components of the resulting graph to see
    // if R intersects a non-trivial component
    // if not, then remove this pair from the list altogether

    std::vector<int32_t> component(num_vertices);
    auto scc = boost::strong_components(H, &component[0]);

    std::vector<std::vector<graph_t::vertex_descriptor>> component_lists(scc);
    boost::build_component_lists(H, scc, &component[0], component_lists);

    for (const auto& list : component_lists) {
      for (const auto& u : list) {
        auto index_u = index[u];
        dbg(OutputType::General, std::print("component({:{}})  =  {:{}}", index_u, width, component[index_u], width));

        // A component with a single vertex is trivial unless the vertex has
        // self-loop.
        if (list.size() == 1) {
          // check for self loop
          auto [e, loop] = boost::edge(u, u, H);

          // self loops are non-trivial strongly connected components
          if (!loop) {
            component[index_u] = TRIVIAL;
          }
        }
      }
    }

    dbg(OutputType::General, std::print("\n"));

    // if a vertex in R lies in a strongly connected component, it may lie
    // on a cycle for some input

    bool trivial = true;
    for (auto i : dynamic_bitset_iterator(pair->right)) {
      dbg(OutputType::General, std::print("component({:{}})  =  {:{}}", i, width, component[i], width));

      if (component[i] != TRIVIAL) {
        trivial = false;
      } else {
        pair->right.reset(i);
      }
    }

    dbg(OutputType::General, std::print("\n"));

    // R does not intersect a non-trivial strongly connected component
    if (trivial) {
      pair = pairs.erase(pair);
    } else {
      ++pair;
    }
  }

  dbg(OutputType::Quiet, Print());
}

// An automaton is universal if the language of the complement is empty.
bool RabinAutomaton::Universal() {
  dbg(OutputType::Quiet, std::print("# Universal()\n\n"));

  const auto width = decimal_digits(num_vertices-1);
  auto universal = false;
  for (auto subset = 0UL; subset < std::exp2(pairs.size()); subset++) {
    graph_t H(graph);
    auto index = boost::get(boost::vertex_index, H);

    // subset \in {0 ... 2^P-1}
    // Specifies a subset of the pairs.
    boost::dynamic_bitset<> pair_set(pairs.size(), subset);
    boost::dynamic_bitset<> clear_set(num_vertices);

    dbg(OutputType::General, std::print("Pair Set: {}\n\n", pair_set));

    // If a vertex appears in any of the specified right conditions, then remove it from the graph.
    const auto pending = ~pair_set;
    for (auto i : dynamic_bitset_iterator(pending)) {
      clear_set = clear_set | pairs[i].right;
    }
    for (auto i : dynamic_bitset_iterator(clear_set)) {
      boost::clear_vertex(boost::vertex(i, H), H);
    }

    auto num_removed = 0UL;
    bool done = false;
    while (!done) {
      done = true;

      for (auto [itr, end] = boost::vertices(H); itr != end; ++itr) {
        if (boost::out_degree(*itr, H) == 0) {
          // Removes all in edges
          // May create additional empty vertices
          boost::clear_vertex(*itr, H);

          // Invalidates iterators so they have to be recreated. Must call
          // clear_vertex() before remove_vertex().
          boost::remove_vertex(*itr, H);

          num_removed++;
          done = false;
          break;
        }
      }
    }

    dbg(OutputType::General, std::print("Clear Set: {}\n\n", clear_set));
    dbg(OutputType::Debug, std::print("Vertices Cleared: {}\n\n", num_removed));

    if (boost::num_vertices(H) == 0) {
      universal = true;
      break;
    }

    // Compute the strongly connected components of the resulting graph to see
    // if a non-trivial strongly connected component intersects every specified
    // left Rabin condition.
    // If so, then the language of the complement automaton is not empty,
    // therefore this automaton cannot be universal.
    std::vector<int32_t> component(boost::num_vertices(H));
    uint32_t scc = boost::strong_components(H, &component[0]);

    std::vector<std::vector<graph_t::vertex_descriptor>> component_lists(scc);
    boost::build_component_lists(H, scc, &component[0], component_lists);

    for (auto& list : component_lists) {
      // A component with a single vertex is trivial unless the vertex has a self-loop.
      if (list.size() == 1) {
        auto& u = list[0];

        // Check for self-loop.
        auto [e, loop] = boost::edge(u, u, H);

        // Self-loops are non-trivial strongly connected components.
        if (!loop) {
          component[index[u]] = TRIVIAL;
        }
      }
    }

    dbg(OutputType::General, {
      std::print("Non-Trivial Components\n");
      for (auto& list : component_lists) {
        for (auto& u : list) {
          auto idx_u = index[u];
          if (component[idx_u] == TRIVIAL) {
            break;
          }

          std::print("component({:{}})  =  {:{}}\n", idx_u, width, component[idx_u], width);
        }
      }

      std::print("\n");
    });

    universal = true;

    // Find any non-trivial component that intersects all specified left Rabin conditions
    for (auto& list : component_lists) {
      auto comp = component[index[list[0]]];
      if (comp == TRIVIAL) {
        continue;
      }

      // If no sets are specified, then any non-trivial component is fine
      if (pair_set.none()) {
        universal = false;
      }

      // Check each rabin pair specified in the pair_set and find a common vertex in the component
      for (auto i : dynamic_bitset_iterator(pair_set)) {
        universal = true;

        for (auto& u : list) {
          auto idx_u = index[u];

          if (pairs[i].left.test(idx_u)) {
            dbg(OutputType::General, {
              std::print("component({:{}})  =  {:{}}\n", idx_u, width, comp, width);
              std::print("{}\n", pairs[i].left);
            });

            universal = false;
            break;
          }
        }

        // The left condition does not intersect the current component
        if (universal) {
          break;
        }
      }

      // This component intersects every specified left Rabin set
      if (!universal) {
        dbg(OutputType::General, std::print("Component {} intersects pair set {}\n", comp, pair_set));

        break;
      }
    }

    if (!universal) {
      dbg(OutputType::General, std::print("The language of the complement automaton is not empty\n\n"));
      break;
    } else {
      dbg(OutputType::General, std::print("No component intersects every left Rabin set\n\n"));
    }
  }

  if (universal) {
    dbg(OutputType::General, std::print("The language of the complement automaton is empty\n"));
  }

  return universal;
}

// Minimization:
// - Remove inaccessible states of M (G is accessible by construction)
// - Compute behavioral equivalence relation for M
// - Merge states with the same behavior
void RabinAutomaton::Minimize() {
  std::vector<boost::dynamic_bitset<>> table;
  table.reserve(num_vertices);
  for (size_type i = 0; i < num_vertices; i++) {
    table.emplace_back(2 * pairs.size());
  }

  // Given a state i, the 2j and 2j+1 bits in the i-th table represent
  // membership in the j-th left and right rabin pair, respectively.

  auto itr = pairs.begin();
  for (auto j = 0UL; j < pairs.size(); j++, itr++) {
    for (auto i = 0UL; i < num_vertices; i++) {
      if (itr->left[i]) {
        table[i].set(2 * j);
      } else if (itr->right[i]) {
        table[i].set(2 * j + 1);
      }
    }
  }

  dbg(OutputType::General, printf("# Minimize()\n\n"));

  const auto width = decimal_digits(num_vertices-1);
  char fmt[MAXLINE];
  std::snprintf(fmt, MAXLINE, "%%%dzd  index = %%%dd  table = ", width, width);

  std::vector<int> equiv(num_vertices);

  // Find identical tables, assign them to the same class, then proceed with
  // normal minimization algorithm.
  std::unordered_map<std::string, size_type> sets;
  for (size_type i = 0; i < num_vertices; i++) {
    std::string s;
    boost::to_string(table[i], s);
    auto itr = sets.find(s);

    if (itr != sets.end()) {
      // This table already exists.
      equiv[i] = itr->second;

      dbg(OutputType::General, std::print("{:{}}  index = {:{}}  table = {}\n", i, width, equiv[i], width, s));
    } else {
      // Store the bitset and make the i-th index the representative.
      sets.insert({s, i});
      equiv[i] = i;

      dbg(OutputType::General, std::print("{:{}}  index = {:{}}  table = {} (*)\n", i, width, equiv[i], width, s));
    }
  }

  dbg(OutputType::General, std::print("\n"));
  std::snprintf(fmt, MAXLINE, "  %%%dd", decimal_digits(num_vertices));

  // std::vector<int> E0(num_vertices);
  std::vector<int> prev(num_vertices);
  std::vector<int> temp(num_vertices);
  std::vector<int> names(num_vertices);

  size_type states = 0;
  size_type rounds = 0;

  auto index = boost::get(boost::vertex_index, graph);
  auto label = boost::get(boost::edge_name, graph);

  bool equivalent = false;
  while (!equivalent) {
    rounds++;

    dbg(OutputType::General, std::print("{:{}}:", "E", binary_digits(num_alphabet-1)));

    // equiv is the initial equivalence class at the start of the round
    for (size_type i = 0; i < num_vertices; i++) {
      temp[i] = -1;
      prev[i] = equiv[i];

      dbg(OutputType::General, std::print("  {:{}}", prev[i], width));
    }

    auto width = binary_digits(num_alphabet-1);
    for (size_type s = 0; s < num_alphabet; s++) {
      dbg(OutputType::General, std::print("\n{:0{}b}:", s, width));

      for (size_type i = 0; i < num_vertices; i++) {
        auto u = boost::vertex(i, graph);

        // Find all edges leaving the current vertex and check if any match
        // the current input symbol.
        for (auto [e_itr, e_end] = boost::out_edges(u, graph);
             e_itr != e_end; ++e_itr) {
          auto v = boost::target(*e_itr, graph);

          if (s != label[*e_itr]) {
            continue;
          }

          // If matching edge is found, update class of vertex to the target.
          auto j = index[v];
          BOOST_ASSERT(j >= 0 && j < num_vertices);

          temp[i] = equiv[j];
          break;
        }

        dbg(OutputType::General, std::print("  {:{}}", temp[i], width));
      }

      dbg(OutputType::General, std::print("\n{:{}}:", "*", binary_digits(num_alphabet-1)));

      size_t current = 0;
      RabinTransitionMap map(num_vertices);

      equivalent = true;
      for (size_type j = 0; j < num_vertices; j++) {
        auto itr = map.find({equiv[j], temp[j]});

        if (itr == map.end()) {
          map.insert({{equiv[j], temp[j]}, j});

          equiv[j] = j;
          names[j] = current++;
        } else {
          equiv[j] = itr->second;
          names[j] = names[equiv[j]];
        }

        temp[j] = -1;
        if (prev[j] != equiv[j]) {
          equivalent = false;
        }

        dbg(OutputType::General, std::print("  {:{}}", equiv[j], width));
      }

      states = map.size();
    }

    dbg(OutputType::General, std::print("\n\n"));
  }

  dbg(OutputType::Quiet, std::print("# Rounds  {}\n", rounds));

  // once behavioral equivalence is computed, merge equivalent states
  if (states == num_vertices) {
    dbg(OutputType::General, std::print("# All states behaviorally unique\n\n"));
    return;
  }

  // create a new graph with only the unique states
  // renumber the states so that the unique states are consecutive
  graph_t H(states);

  label = boost::get(boost::edge_name, H);
  for (auto i = 0UL; i < num_vertices; i++) {
    if (equiv[i] != static_cast<int>(i)) {
      continue;
    }

    auto u = boost::vertex(i, graph);
    for (auto [e_itr, e_end] = boost::out_edges(u, graph);
         e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, graph);
      auto k = index[v];

      assert(k < num_vertices);

      // a is the source and b is the target
      auto a = boost::vertex(names[i], H);
      auto b = boost::vertex(names[k], H);

      auto [e, added] = boost::add_edge(a, b, H);

      // label the transition with the current symbol
      label[e] = boost::get(boost::edge_name, graph, *e_itr);
    }
  }

  dbg(OutputType::General, std::print("# States removed  {}\n\n", num_vertices-states));

  graph = H;

  // make the new rabin pairs
  for (auto& pair : pairs) {
    boost::dynamic_bitset<> L(states);
    boost::dynamic_bitset<> R(states);

    for (std::size_t j = 0; j < num_vertices; j++) {
      if (equiv[j] == static_cast<int>(j)) {
        L.set(names[j], pair.left.test(j));
        R.set(names[j], pair.right.test(j));
      }
    }

    pair.left = L;
    pair.right = R;
  }

  num_vertices = states;
  num_edges = boost::num_edges(graph);

  dbg(OutputType::Quiet, Print());
}

} // namespace omega
