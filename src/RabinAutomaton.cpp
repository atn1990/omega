// RabinAutomaton.cpp
// *  Library to construct and manipulate Rabin automata.

#include "RabinAutomaton.h"
#include "Util.h"

#include <cassert>
#include <cstdio>
#include <iostream>

#include <boost/graph/strong_components.hpp>

namespace omega {

std::ostream& operator<<(std::ostream& os, const RabinPair& p) {
  os << "Left:  " << p.left
     << "\nRight: " << p.right << '\n';
  return os;
}


RabinAutomaton::RabinAutomaton(size_type alphabet, size_type vertices)
  : Automaton(alphabet, vertices) {}

// Build the underlying transition graph from a deterministic transition
// function represented by a map from vertex-symbol pairs to vertices.
void RabinAutomaton::Init(const RabinTransitionMap& map) {
  BOOST_ASSERT(map.size() == num_vertices * num_alphabet);

  // pretty printing of transitions
  char fmt[MAXLINE];
  std::snprintf(fmt, MAXLINE, "(%%%dd,  %%%dd)  -->  %%%dd\n",
      decimal_digits(num_vertices), decimal_digits(num_alphabet), decimal_digits(num_vertices));

  dbg(OutputType::Quiet, printf("# TRANSITIONS\n"));

  auto label = boost::get(boost::edge_name, graph);

  // add edges in sorted order by source vertex and transition symbol
  for (auto i = 0UL; i < num_vertices; i++) {
    for (auto j = 0UL; j < num_alphabet; j++) {
      auto itr = map.find({i, j});

      // the transition system is deterministic
      assert(itr != map.end());
      dbg(OutputType::Quiet, printf(fmt, i, j, itr->second));

      auto [e, b] = boost::add_edge(i, itr->second, graph);
      label[e] = j;
    }
  }

  num_edges = boost::num_edges(graph);
}

void RabinAutomaton::Print() const {
  Automaton::Print();

  printf("# PAIRS\n");
  for (const RabinPair &pair : pairs) {
    std::cout << pair.left << "\n" << pair.right << "\n\n";
  }

  std::flush(std::cout);
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
        dbg(OutputType::Quiet, printf("(%lu,  %u)  ->  %lu\n", v, symbol, t));
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
      dbg(OutputType::Quiet, printf("(%lu,  %u)  ->  %lu (%zd)\n", v, symbol, t, i));
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

  dbg(OutputType::Quiet, printf("# Cycle found (transient = %ld)\n", vertices.size() - V.length()));

  // store cycle in a bitset for comparison against Rabin condition
  boost::dynamic_bitset<> cycle(num_vertices);
  auto index = boost::get(boost::vertex_index, graph);

  // the tail of the vertex list, beginning from start is part of a cycle
  for (auto i = start; i < vertices.size(); i++) {
    cycle.set(index[vertices[i]]);
    dbg(OutputType::Quiet, printf("%lu  ", vertices[i]));
  }

  if (verbose > static_cast<int>(OutputType::Quiet)) {
    std::cout << "\n" << cycle << std::endl;
  }

  // once the cycle stored in a bitset, iterate through Rabin pairs and look
  // for L and R such that L is disjoint from C and R is non-disjoint from C
  for (const auto& pair : pairs) {
    if ((pair.left & cycle).none() && (pair.right & cycle).any()) {
      printf("Input accepted\n");

      if (verbose > static_cast<int>(OutputType::Quiet)) {
        std::cout << "L: " << pair.left << "\nR: " << pair.right << std::endl;
      }

      return;
    }
  }

  // if control reaches here, no accepting Rabin pair was found
  printf("Input rejected\n");
}

/**
  * Removes useless or redundant Rabin pairs. A Rabin pair (L, R) is may be
  * deleted if after removing all vertices from L, there is no non-trivial
  * strongly connected component intersecting R.
 **/
void RabinAutomaton::Clean() {
  dbg(OutputType::Quiet, printf("# Clean\n\n"));

  char fmt[MAXLINE];
  std::snprintf(fmt, MAXLINE, "component(%%%du)  =  %%%dd\n",
      decimal_digits(num_vertices), decimal_digits(num_vertices));

  auto pair = pairs.begin();
  while (pair != pairs.end()) {
    graph_t H(graph);

    // for each vertex in L, remove its corresponding edges
    ITERATE_BITSET(i, pair->left) {
      boost::clear_vertex(boost::vertex(i, H), H);
    }

    if (verbose > static_cast<int>(OutputType::General)) {
      std::cout << pair->left << '\n' << pair->right << "\n\n";
    }

    // Compute the strongly connected components of the resulting graph to see
    // if R intersects a non-trivial component
    // if not, then remove this pair from the list altogether

    std::vector<int32_t> component(num_vertices);
    auto scc = boost::strong_components(H, &component[0]);

    std::vector<std::vector<graph_t::vertex_descriptor>> component_lists(scc);
    boost::build_component_lists(H, scc, &component[0], component_lists);

    for (const auto& list : component_lists) {
      for (const auto& u : list) {
        auto index = boost::get(boost::vertex_index, H, u);
        dbg(OutputType::General, printf(fmt, index, component[index]));

        // A component with a single vertex is trivial unless the vertex has
        // self-loop.
        if (list.size() == 1) {
          // check for self loop
          auto [e, loop] = boost::edge(u, u, H);

          // self loops are non-trivial strongly connected components
          if (!loop) {
            component[index] = TRIVIAL;
          }
        }
      }
    }

    dbg(OutputType::General, printf("\n"));

    // if a vertex in R lies in a strongly connected component, it may lie
    // on a cycle for some input

    bool trivial = true;
    ITERATE_BITSET(i, pair->right) {
      dbg(OutputType::General, printf(fmt, i, component[i]));

      if (component[i] != TRIVIAL) {
        trivial = false;
      } else {
        pair->right.reset(i);
      }
    }

    dbg(OutputType::General, printf("\n"));

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
  char fmt[MAXLINE];
  std::snprintf(fmt, MAXLINE, "component(%%%dzd)  =  %%%dd\n",
           decimal_digits(num_vertices), decimal_digits(num_vertices));

  dbg(OutputType::Quiet, printf("# Universal()\n\n"));

  auto trivial = true;
  for (auto subset = 0UL; subset < std::exp2(pairs.size()); subset++) {
    graph_t H(graph);

    // subset \in 0 -> 2^P - 1
    // Specifies a subset of the pairs.
    boost::dynamic_bitset<> pair_set(pairs.size(), subset);
    boost::dynamic_bitset<> clear_set(num_vertices);

    dbg(OutputType::General, std::cout << "Pair Set: " << pair_set << '\n' << '\n');

    auto i_P = 0U;
    // If a vertex appears in any of the specified right conditions, remove it
    // from the graph.
    for (auto& pair : pairs) {
      if (!pair_set.test(i_P)) {
        clear_set = clear_set | pair.right;
      }

      i_P++;
    }

    std::vector<int32_t> component(num_vertices);
    ITERATE_BITSET(i, clear_set) {
      component[i] = TRIVIAL;
      boost::clear_vertex(boost::vertex(i, H), H);
    }

    // Compute the strongly connected components of the resulting graph to see
    // if a non-trivial strongly connected component intersects every specified
    // left Rabin condition.
    // If so, then the language of the complement automaton is not empty,
    // therefore this automaton cannot be universal.
    uint32_t scc = boost::strong_components(H, &component[0]);

    std::vector<std::vector<graph_t::vertex_descriptor>> component_lists(scc);
    boost::build_component_lists(H, scc, &component[0], component_lists);

    auto index = boost::get(boost::vertex_index, H);

    for (auto& list : component_lists) {
      for (auto& u : list) {
        auto i_U = index[u];
        dbg(OutputType::General, printf(fmt, i_U, component[i_U]));
      }

      // A component with a single vertex is trivial unless vertex has
      // self-loop.
      if (list.size() == 1) {
        auto& u = list[0];
        auto idx = index[u];

        // Check for self-loop.
        auto [e, loop] = boost::edge(u, u, H);

        // Self-loops are non-trivial strongly connected components.
        if (!loop) {
          component[idx] = TRIVIAL;
        }
      }
    }

    if (verbose > static_cast<int>(OutputType::General)) {
      printf("\n");

      for (auto& list : component_lists) {
        for (auto& u : list) {
          auto i_U = index[u];
          printf(fmt, i_U, component[i_U]);
        }
      }

      printf("\n");

      std::cout << "Non-Trivial Components\n";
      for (auto& list : component_lists) {
        for (auto& u : list) {
          auto i_U = index[u];
          if (component[i_U] == TRIVIAL) {
            break;
          }

          printf(fmt, i_U, component[i_U]);
        }
      }

      printf("\n");
    }

    trivial = true;

    // Find any non-trivial component that intersects all specified left Rabin
    // conditions.
    for (auto& list : component_lists) {
      auto c_U = component[index[list[0]]];
      if (c_U == TRIVIAL) {
        continue;
      }

      // If no sets are specified, then any non-trivial component is fine.
      if (pair_set.none()) {
        if (verbose > static_cast<int>(OutputType::General)) {
          printf("Component %d intersects pair set ", c_U);
          std::cout << pair_set << std::endl;
        }

        trivial = false;
        break;
      }

      // check each pair in the set and find a common vertex in the component
      size_t k = 0;
      for (auto& pair : pairs) {
        if (!pair_set.test(k)) {
          k++;
          continue;
        }

        k++;

        trivial = true;
        for (auto& u : list) {
          auto i_U = index[u];

          if (pair.left.test(i_U)) {
            if (verbose > static_cast<int>(OutputType::General)) {
              printf(fmt, i_U, c_U);
              std::cout << pair.left << '\n';
            }

            trivial = false;
            break;
          }
        }

        // The left condition does not intersect the current component.
        if (trivial) {
          break;
        }
      }

      // This component intersects every specified left Rabin set.
      if (!trivial) {
        if (verbose > static_cast<int>(OutputType::General)) {
          printf("Component %d intersects ", c_U);
          std::cout << pair_set << std::endl;
        }

        break;
      }
    }

    if (!trivial) {
      dbg(OutputType::General, printf("The language of the complement automaton is not empty.\n\n"));
      break;
    } else {
      dbg(OutputType::General, printf("No component intersects every left Rabin set.\n\n"));
    }
  }

  if (trivial) {
    dbg(OutputType::General, printf("The language of the complement automaton is empty.\n"));
  }

  return trivial;
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

  char fmt[MAXLINE];
  std::snprintf(fmt, MAXLINE, "%%%dzd  index = %%%dd  table = ",
      decimal_digits(num_vertices), decimal_digits(num_vertices));

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

      dbg(OutputType::General, printf(fmt, i, equiv[i]));
      dbg(OutputType::General, std::cout << table[i] << '\n');
    } else {
      // Store the bitset and make the i-th index the representative.
      sets.insert({s, i});
      equiv[i] = i;

      dbg(OutputType::General, printf(fmt, i, equiv[i]));
      dbg(OutputType::General, std::cout << table[i] << " (*)\n");
    }
  }

  dbg(OutputType::General, printf("\n"));
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

    dbg(OutputType::General, printf("%*s:", binary_digits(num_alphabet), "E"));

    // equiv is the initial equivalence class at the start of the round
    for (size_type i = 0; i < num_vertices; i++) {
      temp[i] = -1;
      prev[i] = equiv[i];

      dbg(OutputType::General, printf(fmt, prev[i]));
    }

    for (size_type s = 0; s < num_alphabet; s++) {
      dbg(OutputType::General,
        printf("\n");
        print_binary(s, binary_digits(num_alphabet));
        printf(":"));

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

        dbg(OutputType::General, printf(fmt, temp[i]));
      }

      dbg(OutputType::General, printf("\n%*s:", binary_digits(num_alphabet), "*"));

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

        dbg(OutputType::General, printf(fmt, equiv[j]));
      }

      states = map.size();
    }

    dbg(OutputType::General, printf("\n\n"));
  }

  dbg(OutputType::Quiet, printf("# Rounds  %llu\n", rounds));

  // once behavioral equivalence is computed, merge equivalent states
  if (states == num_vertices) {
    dbg(OutputType::General, printf("# All states behaviorally unique\n\n"));
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

  dbg(OutputType::General, printf("# States removed  %llu\n\n", num_vertices-states));

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

  if (verbose > static_cast<int>(OutputType::Quiet)) {
    Print();
  }
}

} // namespace omega
