// RabinAutomaton.cpp
// *  Library to construct and manipulate Rabin automata.

#include "RabinAutomaton.h"
#include "Util.h"

#include <cstdio>
#include <functional>
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
// function represented as a dense, row-major table indexed by
// `vertex * num_alphabet + symbol`.
void RabinAutomaton::Init(const std::vector<int_type>& transitions) {
  BOOST_ASSERT(transitions.size() == num_vertices * num_alphabet);

  const auto vertex_width = decimal_digits(num_vertices-1);
  const auto symbol_width = decimal_digits(num_alphabet-1);

  dbg(OutputType::Quiet, std::print("# TRANSITIONS\n"));

  auto label = boost::get(boost::edge_name, graph);

  // add edges in sorted order by source vertex and transition symbol
  for (auto i = 0UL; i < num_vertices; i++) {
    for (auto j = 0UL; j < num_alphabet; j++) {
      auto tgt = transitions[i * num_alphabet + j];
      dbg(OutputType::Quiet, std::print("({:{}}, {:{}})  -->  {:{}}\n", i, vertex_width, j, symbol_width, tgt, vertex_width));

      auto [e, b] = boost::add_edge(i, tgt, graph);
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
//
// The complement of a deterministic Rabin automaton (acceptance: some pair i
// has Inf \cap R_i != {} and Inf \cap L_i == {}) is a Streett condition
// (for every pair i, Inf \cap R_i == {} or Inf \cap L_i != {}). The Rabin
// automaton is universal exactly when this Streett condition has no accepting
// run, which is decided with the classic recursive emptiness test:
//
//   StreettNonEmpty(V):
//     for each non-trivial SCC C of the subgraph induced by V:
//       bad = union of R_i over pairs i with C \cap L_i == {}
//       if bad == {}            -> witness found, return true
//       if StreettNonEmpty(C \ bad) -> return true
//     return false
//
// The computation is performed over the original vertex identifiers (using a
// filtered view of the graph) so that the indices always line up with the
// bitsets stored in the Rabin pairs.
bool RabinAutomaton::Universal() {
  dbg(OutputType::General, std::print("# Universal()\n\n"));

  auto index = boost::get(boost::vertex_index, graph);

  const int_type n = num_vertices;

  // Build a compressed sparse row (CSR) adjacency once so the emptiness test
  // below can explore the graph without re-iterating Boost edge lists at every
  // recursion level. out_begin[v]..out_begin[v+1] indexes v's successors in
  // `targets`; self_loop records vertices with an edge to themselves.
  std::vector<int_type> out_begin(n + 1, 0);
  for (int_type v = 0; v < n; v++) {
    out_begin[v + 1] =
      out_begin[v] + boost::out_degree(boost::vertex(v, graph), graph);
  }

  std::vector<int_type> targets(out_begin[n]);
  boost::dynamic_bitset<> self_loop(n);
  {
    std::vector<int_type> pos(out_begin.begin(), out_begin.end() - 1);
    for (int_type v = 0; v < n; v++) {
      for (auto [itr, end] = boost::out_edges(boost::vertex(v, graph), graph);
           itr != end; ++itr) {
        auto w = index[boost::target(*itr, graph)];
        targets[pos[v]++] = w;
        if (w == v) {
          self_loop.set(v);
        }
      }
    }
  }

  // Scratch state shared across all recursion levels. Tarjan visitation is
  // marked with a monotonically increasing epoch so the O(n) arrays never need
  // to be cleared between calls; `onstack` is always left empty when a call
  // returns. This keeps each level proportional to the active sub-automaton
  // rather than the whole graph.
  std::vector<int_type> ord(n);
  std::vector<int_type> low(n);
  std::vector<int_type> stamp(n, 0);
  boost::dynamic_bitset<> onstack(n);
  std::vector<int_type> dfs_stack;
  std::vector<int_type> ei_stack;
  std::vector<int_type> scc_stack;
  int_type epoch = 0;

  // Iterative Tarjan restricted to the active vertices and the edges between
  // them. Collects the strongly connected components into `comps`.
  auto strong_components_active =
    [&](const boost::dynamic_bitset<>& active,
        std::vector<std::vector<int_type>>& comps) {
    comps.clear();
    ++epoch;
    int_type counter = 0;

    for (auto root : dynamic_bitset_iterator(active)) {
      if (stamp[root] == epoch) {
        continue;
      }

      dfs_stack.clear();
      ei_stack.clear();
      dfs_stack.push_back(root);
      ei_stack.push_back(out_begin[root]);
      stamp[root] = epoch;
      ord[root] = low[root] = counter++;
      scc_stack.push_back(root);
      onstack.set(root);

      while (!dfs_stack.empty()) {
        auto v = dfs_stack.back();
        auto ei = ei_stack.back();
        bool descended = false;

        while (ei < out_begin[v + 1]) {
          auto w = targets[ei++];
          if (!active[w]) {
            continue;
          }
          if (stamp[w] != epoch) {
            ei_stack.back() = ei;
            stamp[w] = epoch;
            ord[w] = low[w] = counter++;
            scc_stack.push_back(w);
            onstack.set(w);
            dfs_stack.push_back(w);
            ei_stack.push_back(out_begin[w]);
            descended = true;
            break;
          } else if (onstack[w] && ord[w] < low[v]) {
            low[v] = ord[w];
          }
        }

        if (descended) {
          continue;
        }

        // v is fully explored: close its SCC if it is a root, then backtrack.
        if (low[v] == ord[v]) {
          std::vector<int_type> component;
          int_type w;
          do {
            w = scc_stack.back();
            scc_stack.pop_back();
            onstack.reset(w);
            component.push_back(w);
          } while (w != v);
          comps.push_back(std::move(component));
        }

        dfs_stack.pop_back();
        ei_stack.pop_back();
        if (!dfs_stack.empty()) {
          auto p = dfs_stack.back();
          if (low[v] < low[p]) {
            low[p] = low[v];
          }
        }
      }
    }
  };

  // Recursive Streett emptiness test on the set of active vertices.
  std::function<bool(const boost::dynamic_bitset<>&)> non_empty =
    [&](const boost::dynamic_bitset<>& active) -> bool {
    if (active.none()) {
      return false;
    }

    std::vector<std::vector<int_type>> comps;
    strong_components_active(active, comps);

    // Reused across components to avoid repeated full-width allocations.
    boost::dynamic_bitset<> C(n);

    for (auto& member : comps) {
      // A single-vertex component is non-trivial only if it has a self-loop.
      if (member.size() == 1 && !self_loop[member.front()]) {
        continue;
      }

      C.reset();
      for (auto v : member) {
        C.set(v);
      }

      // Pairs whose left condition is disjoint from the component cannot be
      // satisfied here, so their right-condition states must be avoided.
      boost::dynamic_bitset<> bad(n);
      for (const auto& pair : pairs) {
        if (!C.intersects(pair.left)) {
          bad |= (pair.right & C);
        }
      }

      // No states need to be removed: the whole component is a witness for a
      // non-empty complement language.
      if (bad.none()) {
        return true;
      }

      if (non_empty(C & ~bad)) {
        return true;
      }
    }

    return false;
  };

  // Restrict the search to the part of the automaton reachable from the
  // initial state (vertex 0).
  boost::dynamic_bitset<> reachable(num_vertices);
  if (num_vertices > 0) {
    std::vector<graph_t::vertex_descriptor> stack = {boost::vertex(0, graph)};
    reachable.set(0);
    while (!stack.empty()) {
      auto u = stack.back();
      stack.pop_back();
      for (auto [itr, end] = boost::out_edges(u, graph); itr != end; ++itr) {
        auto v = boost::target(*itr, graph);
        if (!reachable[index[v]]) {
          reachable.set(index[v]);
          stack.push_back(v);
        }
      }
    }
  }

  auto universal = !non_empty(reachable);

  if (universal) {
    dbg(OutputType::General, std::print("The language of the complement automaton is empty\n"));
  } else {
    dbg(OutputType::General, std::print("The language of the complement automaton is not empty\n"));
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
  [[maybe_unused]] size_type rounds = 0;

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
