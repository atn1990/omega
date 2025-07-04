// BuchiAutomaton.cpp
// * Implementation of the BuchiAutomaton class.

#include "BuchiAutomaton.h"
#include "RabinAutomaton.h"
#include "Safra.h"
#include "Util.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <queue>
#include <source_location>
#include <unordered_set>
#include <unordered_map>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/visitors.hpp>

// TODO: See if any macros from here are useful?
// #include <boost/graph/graph_utility.hpp>
// #include <boost/graph/iteration_macros.hpp>

namespace logger {
  bool enabled{true};
} // namespace logger

namespace omega {

BuchiAutomaton::BuchiAutomaton(size_type num_alphabet, size_type num_vertices)
  : Automaton(num_alphabet, num_vertices), initial_states(num_vertices), final_states(num_vertices) {}

/**
  * Given a non-deterministic transition function, fill in the initialized
  * graph with the corresponding edges.
  *
  * map is an unordered_multimap representing the edges
  * (vertex * symbol -> vertex) of the non-deterministic Büchi automaton.
  *
  * Assumes the underlying graph of the Büchi automaton is initialized with the
  * correct number of states.
 **/
void BuchiAutomaton::Init(const TransitionMap& map) {
  auto label = boost::get(boost::edge_name, graph);
  auto state = boost::get(boost::vertex_name, graph);

  for (auto i = 0U; i < num_vertices; i++) {
    auto u = boost::vertex(i, graph);

    for (auto j = 0U; j < num_alphabet; j++) {
      // find all edges leaving vertex i under symbol j
      for (auto [itr, end] = map.equal_range({i, j}); itr != end; ++itr) {
        // the second element of the iterator is the new state
        auto v = boost::vertex(itr->second, graph);

        auto [e, added] = boost::add_edge(u, v, graph);
        label[e] = j;
      }
    }

    if (initial_states[i] && final_states[i]) {
      state[u] = NodeType::Both;
    } else if (initial_states[i]) {
      state[u] = NodeType::Initial;
    } else if (final_states[i]) {
      state[u] = NodeType::Final;
    } else {
      state[u] = NodeType::None;
    }
  }

  num_edges = boost::num_edges(graph);
}

/**
  * Pretty printing of Büchi automata.
  *
  * Supports printing to a file format (buchi.aut) for storage.
 **/
void BuchiAutomaton::Print() const {
  char fmt[MAXLINE];
  snprintf(fmt, MAXLINE, "  %%%dzd", decimal_digits(num_vertices));

  if (verbose == static_cast<int>(OutputType::Outfile)) {
    printf("# BUCHI\n");
    printf("# RABIN SIZE: ?\n");
    printf("# RABIN TRANSITIONS: ?\n");

    printf("# STATES\n%llu\n", num_vertices);
    printf("# ALPHABET\n%llu\n", num_alphabet);
    printf("# TRANSITIONS\n%llu\n", num_edges);

    printf("# BEGIN TRANSITIONS\n");
    Automaton::Print();
    printf("# END TRANSITIONS\n");

    printf("# INITIAL\n");
    std::cout << initial_states << "\n";

    printf("# FINAL\n");
    std::cout << final_states << "\n";

    printf("# BUCHI EOF");
  } else {
    printf("# STATES      %llu\n", num_vertices);
    printf("# ALPHABET    %llu\n", num_alphabet);
    printf("# TRANSITIONS %llu\n\n", num_edges);

    Automaton::Print();

    printf("# INITIAL\n");
    if (verbose > static_cast<int>(OutputType::General)) {
      ITERATE_BITSET(i, initial_states) {
        printf(fmt, i);
      }

      if (initial_states.none()) {
        printf("Empty");
      }

      std::cout << "\n";
    } else {
      std::cout << initial_states << "\n";
    }

    printf("# FINAL\n");
    if (verbose > static_cast<int>(OutputType::General)) {
      ITERATE_BITSET(i, final_states) {
        printf(fmt, i);
      }

      if (final_states.none()) {
        printf("Empty");
      }

      std::cout << "\n";
    } else {
      std::cout << final_states << "\n";
    }
  }

  std::cout << std::endl << std::endl;
}

// Helper function to update the internal variables of the Büchi automaton
// after the underlying graph has been modified to ensure algorithms still
// work correctly.
//
// Uses the vertex_name property of the vertices to correctly keep track of
// the initial and final states even if their id has changed.
void BuchiAutomaton::Resize() {
  Automaton::Resize();

  initial_states.resize(num_vertices);
  initial_states.reset();

  final_states.resize(num_vertices);
  final_states.reset();

  auto state = boost::get(boost::vertex_name, graph);

  auto i = 0U;
  for (auto [itr, end] = boost::vertices(graph); itr != end; ++itr) {
    if (state[*itr] == NodeType::Both) {
      initial_states.set(i);
      final_states.set(i);
    } else if (state[*itr] == NodeType::Initial) {
      // TODO: Are initial states also final?
      initial_states.set(i);
    } else if (state[*itr] == NodeType::Final) {
      final_states.set(i);
    }

    i++;
  }
}

// Removes all useless vertices from the underlying transition graph.
// A vertex is useless if it has no outgoing edges, or if its outgoing edges
// are only to useless vertices.
void BuchiAutomaton::Clean() {
  auto done = false;
  auto num_removed = 0U;

  while (!done) {
    done = true;

    for (auto [itr, end] = boost::vertices(graph); itr != end; ++itr) {
      if (boost::out_degree(*itr, graph) != 0) {
        continue;
      }

      // Removes all in edges. May create additional empty vertices.
      boost::clear_vertex(*itr, graph);

      // Invalidates iterators so they have to be recreated. Must call
      // clear_vertex() before remove_vertex().
      boost::remove_vertex(*itr, graph);

      num_removed++;
      done = false;
      break;
    }
  }

  dbg(OutputType::Debug, printf("Vertices Removed: %u\n\n", num_removed));

  // The graph may have changed so update the state of the automaton.
  Resize();
}

// custom dfs visitor to compute the accessible part of an automaton
class dfs_reachable_visitor : public boost::default_dfs_visitor {
public:
  dfs_reachable_visitor(boost::dynamic_bitset<>& set) : reachable(set) {}

  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, const Graph &) {
    // dbg(OutputType::Debug, std::cout << "discovered: " << u << "\n");
    reachable.set(u);
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, const Graph&) const {
    // dbg(OutputType::Debug, std::cout << "finished: " << u << "\n\n");
  }

  template <typename Edge, typename Graph>
  void examine_edge(Edge e, const Graph&) const {
    // dbg(OutputType::Debug, std::cout << "examined: " << e << "\n");
  }

  boost::dynamic_bitset<>& reachable;
};

// Identify and remove inaccessible states from the underlying graph.
void BuchiAutomaton::Reachable() {
  // visitor keeps track of all discovered vertices
  boost::dynamic_bitset<> reachable(num_vertices);
  dfs_reachable_visitor visitor(reachable);

  auto index = boost::get(boost::vertex_index, graph);
  // only search the connected components of initial states
  ITERATE_BITSET(i, initial_states) {
    auto v = boost::vertex(i, graph);

    color_map_t color_map(num_vertices, index);
    boost::depth_first_visit(graph, v, visitor, color_map);
  }

  if (verbose > static_cast<int>(OutputType::Debug)) {
    std::cout << "Reachable: " << reachable << "\n\n";
  }

  reachable.flip();

  if (reachable.any()) {
    ITERATE_BITSET(i, reachable) {
      auto v = boost::vertex(i, graph);
      boost::clear_vertex(v, graph);
    }

    Clean();
  }

  if (verbose > static_cast<int>(OutputType::Debug)) {
    Print();
  }
}


// Project away the most significant symbol from the edges of the underlying
// transition graph.
// A label is a sequence of bits corresponding to input on a set of tracks.
// Removing a track may result in a non-deterministic transition system.
void BuchiAutomaton::ProjectLabel() {
  if (verbose > static_cast<int>(OutputType::General)) {
    printf("# ProjectLabel()\n");
  }

  // TODO: Should num_alphabet be a power of two?
  BOOST_ASSERT(num_alphabet % 2 == 0);

  // The alphabet of the automaton is cut in half.
  num_alphabet >>= 1;

  // Create a mask to erase the most significant bit of the label on an edge.
  // Note: only works if num_alphabet is a power of two.
  size_type mask = num_alphabet - 1;

  auto index = boost::get(boost::vertex_index, graph);
  auto label = boost::get(boost::edge_name, graph);

  for (auto [v_itr, v_end] = boost::vertices(graph); v_itr != v_end; ++v_itr) {
    // Explore all outgoing edges of the current vertex.
    auto done = false;
    while (!done) {
      done = true;

      std::unordered_set<int_pair> edges;
      for (auto [e_itr, e_end] = boost::out_edges(*v_itr, graph); e_itr != e_end; ++e_itr) {
        auto t = boost::target(*e_itr, graph);
        auto symbol = label[*e_itr] & mask;
        label[*e_itr] = symbol;

        auto itr = edges.find({index[t], symbol});

        // Update the label on the edge and record it to prevent duplicates
        if (itr == edges.end()) {
          edges.insert({index[t], symbol});
        } else {
          // Remove the duplicated edge.
          boost::remove_edge(*e_itr, graph);
          done = false;
          break;
        }
      }
    }
  }

  // by removing duplicated edges, no vertices become inaccessible
  num_edges = boost::num_edges(graph);

  if (verbose > static_cast<int>(OutputType::General)) {
    Print();
  }
}

// Generate a representation of the underlying transition graph as an
// unordered hash multimap.
//
// The keys are vertex-symbol pairs and the values are vertices.
std::unique_ptr<TransitionMap> BuchiAutomaton::Map() const {
  auto map = std::make_unique<TransitionMap>();

  auto index = boost::get(boost::vertex_index, graph);
  auto symbol = boost::get(boost::edge_name, graph);

  for (auto [itr, end] = boost::edges(graph); itr != end; ++itr) {
    auto u = boost::source(*itr, graph);
    auto v = boost::target(*itr, graph);

    map->insert({{index[u], symbol[*itr]}, index[v]});
  }

  return map;
}

// Display a path in a graph ending at a target vertex v using the vertices in
// the predecessor list p, generated from a breadth first search traversal of
// the underlying transition graph.
void BuchiAutomaton::PrintPath(
    std::vector<graph_t::vertex_descriptor>& pred,
    graph_t::vertex_descriptor& v,
    std::vector<uint8_t>& path,
    char* fmt_s,
    char* fmt_t) {
  auto index  = boost::get(boost::vertex_index, graph);
  auto symbol = boost::get(boost::edge_name, graph);

  if (pred[v] == v) {
    return;
  }

  PrintPath(pred, pred[v], path, fmt_s, fmt_t);

  for (auto [itr, end] = boost::in_edges(v, graph); itr != end; ++itr) {
    auto s = boost::source(*itr, graph);
    if (index[s] != index[pred[v]]) {
      continue;
    }

    path.push_back(symbol[*itr]);

    printf(fmt_s, index[s]);
    print_binary(symbol[*itr], binary_digits(num_alphabet));
    printf(fmt_t, index[v]);

    return;
  }
}

// Find a cycle within a strongly connected component containing the source
// vertex s.
void BuchiAutomaton::FindCycle(
    graph_t::vertex_descriptor& s,
    std::vector<int>& component,
    std::vector<uint8_t>& path,
    char* fmt_s,
    char* fmt_t) {
  auto vertex_index = boost::get(boost::vertex_index, graph);
  auto symbol = boost::get(boost::edge_name, graph);

  auto [e, loop] = boost::edge(s, s, graph);
  if (loop) {
    path.push_back(symbol[e]);

    printf(fmt_s, vertex_index[s]);
    print_binary(symbol[e], binary_digits(num_alphabet));
    printf(fmt_t, vertex_index[s]);

    return;
  }

  // Since s is in a strongly connected component and it doesn't have a
  // self-loop, there must exist some neighbor t in the same component.
  graph_t::vertex_descriptor t;

  for (auto [itr, end] = boost::out_edges(s, graph); itr != end; ++itr) {
    t = boost::target(*itr, graph);

    if (component[vertex_index[t]] == component[vertex_index[s]]) {
      path.push_back(symbol[*itr]);

      printf(fmt_s, vertex_index[s]);
      print_binary(symbol[*itr], binary_digits(num_alphabet));
      printf(fmt_t, vertex_index[t]);

      break;
    }
  }

  // Predecessor Map
  std::vector<graph_t::vertex_descriptor> pred(num_vertices);

  auto visitor =
    boost::make_bfs_visitor(
        boost::record_predecessors(&pred[0], boost::on_tree_edge()));

  boost::breadth_first_search(graph, t, boost::visitor(visitor));

  PrintPath(pred, s, path, fmt_s, fmt_t);
}

// Checks if the language recognized by the automaton is empty.
bool BuchiAutomaton::Empty() {
  dbg(OutputType::Quiet, printf("# EMPTY\n\n"));
  if (final_states.none()) {
    // If there are no final states, the language is empty.
    dbg(OutputType::Quiet, printf("The language is empty\n"));
    return true;
  }

  // Generate the strongly connected components using Tarjan's algorithm.
  std::vector<int> component(num_vertices);
  auto scc = boost::strong_components(graph, &component[0]);
  dbg(OutputType::Debug, printf("SCC: %d\n", scc));

  // Generate the component lists.
  std::vector<std::vector<graph_t::vertex_descriptor>> component_lists(scc);
  boost::build_component_lists(graph, scc, &component[0], component_lists);

  char fmt[MAXLINE];
  snprintf(fmt, MAXLINE, "component(%%%du)  =  %%%dd\n",
      decimal_digits(num_vertices), decimal_digits(scc));

  auto vertex_index = boost::get(boost::vertex_index, graph);

  // Identify the non-trivial components.
  for (auto& list : component_lists) {
    bool trivial = true;

    for (auto& u : list) {
      auto k = vertex_index[u];
      dbg(OutputType::Debug, printf(fmt, k, component[k]));

      // A component is non-trivial if it contains a final state.
      if (final_states[k]) {
        trivial = false;
      }
    }

    if (trivial) {
      // Mark all vertices in this component as trivial.
      for (auto& u : list) {
        component[vertex_index[u]] = TRIVIAL;
      }
    } else if (list.size() == 1) {
      auto& u = list[0];
      auto [e, loop] = boost::edge(u, u, graph);

      // A component of exactly one final state is trivial unless the state has
      // a self loop.
      if (!loop) {
        component[vertex_index[u]] = TRIVIAL;
      }
    }
  }

  if (verbose > static_cast<int>(OutputType::General)) {
    std::cout << "\nNon-Trivial Components\n";

    bool all_trivial = true;

    // Print the non-trivial components.
    for (auto& list : component_lists) {
      for (auto& u : list) {
        auto k = vertex_index[u];
        if (component[k] != TRIVIAL) {
          all_trivial = false;
          printf(fmt, k, component[k]);
        }
      }
    }

    if (all_trivial) {
      std::cout << "All components are trivial\n";
    }

    std::cout << '\n';
  }

  char fmt_s[MAXLINE];
  snprintf(fmt_s, MAXLINE, "(%%%du, ", decimal_digits(num_vertices));

  char fmt_t[MAXLINE];
  snprintf(fmt_t, MAXLINE, ")  -->  %%%du\n", decimal_digits(num_vertices));

  bool trivial = true;
  ITERATE_BITSET(i, initial_states) {
    auto u = boost::vertex(i, graph);

    // Predecessor Map
    std::vector<graph_t::vertex_descriptor> pred(num_vertices);
    color_map_t color_map(num_vertices, vertex_index);

    auto visitor =
      boost::make_bfs_visitor(
          boost::record_predecessors(&pred[0], boost::on_tree_edge()));

    // Do a BFS from each of the initial vertices.
    boost::breadth_first_search(
        graph, u, boost::visitor(visitor).color_map(color_map));

    ITERATE_BITSET(j, final_states) {
      auto v = boost::vertex(j, graph);

      if (component[j] == TRIVIAL) {
        continue;
      }

      // If color[v] == white, then v was unvisited.
      if (boost::get(color_map, v) == color_t::white()) {
        continue;
      }

      // There is a path from an initial state to a final state in a
      // non-trivial strongly connected component.
      if (verbose > static_cast<int>(OutputType::General)) {
        std::vector<uint8_t> path;

        // Output the path from the initial state to the final state.
        printf("%zd  -->  %zd\n", i, j);
        PrintPath(pred, v, path, fmt_s, fmt_t);

        // Output a cycle in the strongly connected component.
        printf("\n%zd  -->  %zd\n", j, j);
        FindCycle(v, component, path, fmt_s, fmt_t);
        printf("\n");

        auto tracks = binary_digits(num_alphabet);
        // Output the labels of the edges visited.
        for (auto k = 0; k < tracks; k++) {
          // Show the phantom zero.
          if (k < tracks - 1) {
            printf("0  ");
          } else {
            printf("   ");
          }

          for (auto l : path) {
            printf("%llu  ", MAP(l, k));
          }

          printf("\n");
        }

        printf("\n");
      }

      trivial = false;
    }

    if (!trivial) {
      break;
    }
  }

  if (verbose > static_cast<int>(OutputType::General)) {
    std::cout << std::endl;
  }

  return trivial;
}

inline void who(const std::source_location sloc = std::source_location::current()) {
  if (logger::enabled) {
    std::cout << sloc.function_name() << '\n';
  }
}

struct SafraTreeHash {
  size_t operator()(const std::shared_ptr<SafraTree>& tree) const {
    // who();
    return tree->hash_value();
  }
};

struct SafraTreeEq {
  bool operator()(const std::shared_ptr<SafraTree>& t1,
                  const std::shared_ptr<SafraTree>& t2) const {
    // who();
    if (t1 != nullptr && t2 != nullptr) {
      return (*t1) == (*t2);
    } else {
      return t1 == t2;
    }
  }
};

template <class Name>
class label_writer {
  public:
    label_writer(Name _name) : name(_name) {}

    template <class VertexOrEdge>
    void operator()(std::ostream& os, const VertexOrEdge& v) const {
      os << "[label=\"" << static_cast<int_type>(name[v]) << "\"]";
    }

  private:
    Name name;
};

// Converts a non-deterministic Büchi automaton into a deterministic Rabin
// automaton using Safra's determinization algorithm.

// Safra's algorithm consists of the following six steps:
//   - unmark
//   - update
//   - create
//   - horizontal_merge
//   - kill_empty
//   - vertical_merge

// Starting from the initial tree, these steps are iterated on each tree for
// each input symbol until no new tree is produced.

// A queue of SafraTree pointers keeps track of the next tree to compute:
//   - pop a tree from the front of the queue and run the six steps of the
//     algorithm for every input symbol
//   - if a tree is unique, hash it, give it a unique index number, and push
//     it to the end of the queue

// A std::unordered_multimap keeps track of the transitions for the input
// automaton, as well as the computed Rabin automaton.

// Every time a non-empty tree is found, check if it is in the hash table:
//   - if it is, retrieve the stored tree, and store a new transition from the
//     index of the old tree and the current input symbol to the index of the
//     stored tree
//   - otherwise, use the lowest index number available and store a new
//     transition as before

// Once the algorithm terminates, convert the hashed transitions into a graph
// with the hashed trees as vertices. Extract the Rabin pairs from the stored
// trees.
std::unique_ptr<RabinAutomaton> BuchiAutomaton::Determinize(const TransitionMap &map, DeterminizeOpts opts) const {
  std::unordered_set<std::shared_ptr<SafraTree>, SafraTreeHash, SafraTreeEq> set;
  std::queue<std::shared_ptr<SafraTree>> queue;

  RabinTransitionMap rabin_map;

  auto empty_index = -1L;

  auto num_empty = 0UL;
  auto num_trees = 1UL;

  auto initial = std::make_shared<SafraTree>(num_vertices, initial_states, final_states);

  // The initial tree has index 0
  if (verbose > static_cast<int>(OutputType::General)) {
    printf("# Determinize()\n");
    printf("Index:  %u\n", initial->index);
    initial->PrintTree();
  }

  // hash the initial tree and push a copy onto the queue
  set.insert(initial);
  queue.push(std::make_shared<SafraTree>(*initial));

  // Guards the queue and condition variable.
  // std::mutex queue_mutex;
  // std::condition_variable cv;

  // Guards access to set, list, and rabin_map.
  // std::mutex data_mutex;

  // auto ready = true;
  // auto num_waiting = 0UL;

  // auto f =
  //   [&]() -> void {
  //     while (true) {
  //       std::unique_lock<std::mutex> queue_lock(queue_mutex);

  //       num_waiting++;
  //       if (num_waiting == num_threads && queue.empty()) {
  //         ready = true;
  //         queue_lock.unlock();
  //         cv.notify_one();
  //         return;
  //       }

  //       cv.wait(queue_lock, [&ready]() { return ready; });

  //       // queue_mutex is locked and ready == true
  //       // If queue is empty, then no more work to do.
  //       if (queue.empty()) {
  //         queue_lock.unlock();
  //         cv.notify_one();
  //         return;
  //       }

  //       auto T = queue.front();
  //       queue.pop();

  //       num_waiting--;
  //       if (!queue.empty()) {
  //         ready = true;
  //         queue_lock.unlock();
  //         cv.notify_one();
  //       } else {
  //         ready = false;
  //         queue_lock.unlock();
  //       }

  //       {
  //         std::ostringstream buf;
  //         buf << "Thread " << std::this_thread::get_id()
  //             << " processing tree " << T->index << std::endl;
  //         std::cerr << buf.str();
  //       }

  //       if (verbose > OutputType::General) {
  //         std::ostringstream buf;
  //         buf << "Index:  " << T->index << '\n';
  //         T->PrintTree(buf);
  //         std::cout << buf.str();
  //       }

  //       // save time by unmarking the tree once
  //       T->Unmark();
  //       dbg(OutputType::Debug,
  //           std::ostringstream buf;
  //           buf << "# Unmark()\n";
  //           T->PrintTree(buf);
  //           std::cout << buf.str());

  //       // If Update and Create are swapped, children are created only once
  //       // per new tree.
  //       if (opts.SwapUpdateCreate) {
  //         T->Create(final_states);
  //         dbg(OutputType::Debug,
  //             std::ostringstream buf;
  //             buf << "# Create()\n";
  //             T->PrintTree(buf);
  //             std::cout << buf.str());
  //       }

  //       for (auto i = 0U; i < num_alphabet; i++) {
  //         auto U = std::make_shared<SafraTree>(*T);
  //         num_trees++;

  //         U->Update(map, i);
  //         dbg(OutputType::Debug,
  //             std::ostringstream buf;
  //             buf << "# Update() f(_, ";
  //             print_binary(i, binary_digits(num_alphabet), buf);
  //             buf << ")\n";
  //             U->PrintTree(buf);
  //             std::cout << buf.str());

  //         // The update and create steps are interchangable.
  //         if (!opts.SwapUpdateCreate) {
  //           U->Create(final_states);
  //           dbg(OutputType::Debug,
  //               std::ostringstream buf;
  //               buf << "# Create()\n";
  //               U->PrintTree(buf);
  //               std::cout << buf.str());
  //         }

  //         U->HorizontalMerge();
  //         dbg(OutputType::Debug,
  //             std::ostringstream buf;
  //             buf << "# HorizontalMerge()\n";
  //             U->PrintTree(buf);
  //             std::cout << buf.str());

  //         U->KillEmpty();
  //         dbg(OutputType::Debug,
  //             std::ostringstream buf;
  //             buf << "# KillEmpty()\n";
  //             U->PrintTree(buf);
  //             std::cout << buf.str());

  //         // Delete empty trees.
  //         if (U->root == nullptr) {
  //           std::lock_guard<std::mutex> lock(data_mutex);
  //           // data_mutex.lock();
  //           // mutex.lock();
  //           // rw.lock();

  //           if (num_empty == 0) {
  //             // boost::dynamic_bitset<> L(num_vertices);

  //             // U->names.reset();
  //             // U->marked.reset();

  //             // Create a single root node with empty label and unmarked.
  //             // U->Init(L, ~L);

  //             U->index = list.size();
  //             empty_index = U->index;

  //             set.insert(U);
  //             list.push_back(U);

  //             // The empty tree is an absorbing state with a self-loop under
  //             // all transition symbols
  //             for (auto j = 0UL; j < num_alphabet; j++) {
  //               rabin_map.insert({{empty_index, j}, empty_index});
  //             }

  //             rabin_map.insert({{T->index, i}, empty_index});
  //             num_empty = num_alphabet + 1;
  //             // data_mutex.unlock();

  //             dbg(OutputType::General, printf("Unique:  %u\n\n", U->index));
  //           } else {
  //             rabin_map.insert({{T->index, i}, empty_index});
  //             num_empty++;
  //             // data_mutex.unlock();

  //             dbg(OutputType::General, printf("\nEmpty:  %llu\n\n", empty_index));
  //           }

  //           continue;
  //         }

  //         U->VerticalMerge();
  //         dbg(OutputType::Debug,
  //             std::ostringstream buf;
  //             buf << "# VerticalMerge()\n";
  //             U->PrintTree(buf);
  //             std::cout << buf.str());

  //         data_mutex.lock();
  //         auto itr = set.find(U);

  //         if (itr == set.end()) {
  //           // A new tree
  //           // - Add it to the unordered_set of trees
  //           // - Push it on the queue of pending trees
  //           // - Add the corresponding transition
  //           U->index = list.size();

  //           list.push_back(U);
  //           set.insert(U);

  //           rabin_map.insert({{T->index, i}, U->index});
  //           data_mutex.unlock();

  //           queue_mutex.lock();
  //           queue.push(std::make_shared<SafraTree>(*U));

  //           ready = true;
  //           queue_mutex.unlock();
  //           cv.notify_one();

  //           dbg(OutputType::General, printf("Unique:  %u\n\n", U->index));
  //         } else {
  //           // same as a previously computed tree
  //           // add the corresponding transition
  //           rabin_map.insert({{T->index, i}, (*itr)->index});
  //           auto index = (*itr)->index;
  //           data_mutex.unlock();

  //           dbg(OutputType::General, printf("Index:  %u\n\n", index));
  //         }
  //       }
  //     }
  //   };

  // std::vector<std::thread> threads(num_threads);
  // for (auto i = 0U; i < num_threads; i++) {
  //   threads[i] = std::thread(f);
  // }

  // for (auto i = 0U; i < num_threads; i++) {
  //   threads[i].join();
  // }

  while (!queue.empty()) {
    auto T = queue.front();
    queue.pop();

    if (verbose > static_cast<int>(OutputType::General)) {
      printf("Index:  %u\n", T->index);
      T->PrintTree();
    }

    // save time by unmarking the tree once
    T->Unmark();
    if (verbose > static_cast<int>(OutputType::General)) {
      printf("# Unmark()\n");
      T->PrintTree();
    }

    // if update and create are swapped, children are created once per new tree
    if (opts.SwapUpdateCreate) {
      T->Create(final_states);
      if (verbose > static_cast<int>(OutputType::General)) {
        printf("# Create()\n");
        T->PrintTree();
      }
    }

    for (auto i = 0U; i < num_alphabet; i++) {
      auto U = std::make_shared<SafraTree>(*T);
      num_trees++;

      U->Update(map, i);
      if (verbose > static_cast<int>(OutputType::General)) {
        printf("# Update() f(_, ");
        print_binary(i, binary_digits(num_alphabet));
        printf(")\n");
        U->PrintTree();
      }

      // update and create steps are interchangable
      if (!opts.SwapUpdateCreate) {
        U->Create(final_states);
        if (verbose > static_cast<int>(OutputType::General)) {
          printf("# Create()\n");
          U->PrintTree();
        }
      }

      U->HorizontalMerge();
      if (verbose > static_cast<int>(OutputType::General)) {
        printf("# HorizontalMerge()\n");
        U->PrintTree();
      }

      U->KillEmpty();
      if (verbose > static_cast<int>(OutputType::General)) {
        printf("# KillEmpty()\n");
        U->PrintTree();
      }

      // delete empty trees
      if (U->root == nullptr) {
        if (num_empty == 0) {
          U->names.reset();
          // U->marked.reset();

          // create a single root node with empty label and unmarked
          // boost::dynamic_bitset<> L(num_vertices);
          // U->Init(L, ~L);

          empty_index = set.size();
          U->index = empty_index;

          set.insert(U);

          // The empty tree is an absorbing state with a self-loop under all
          // transition symbols.
          for (auto j = 0U; j < num_alphabet; j++) {
            rabin_map.insert({{empty_index, j}, (int_type) empty_index});
          }

          num_empty = num_alphabet;
          num_trees += num_alphabet;

          dbg(OutputType::General, printf("# Unique  %u\n\n", U->index));
        } else {
          dbg(OutputType::General, printf("# Empty  %ld\n\n", empty_index));
        }

        rabin_map.insert({{T->index, i}, (int_type) empty_index});
        num_empty++;

        continue;
      }

      U->VerticalMerge();
      if (verbose > static_cast<int>(OutputType::General)) {
        printf("# VerticalMerge()\n");
        U->PrintTree();
      }

      auto itr = set.find(U);
      if (itr == set.end()) {
        // unique tree
        //  - add it to the hash table of trees
        //  - push it to the queue of pending trees
        //  - add the corresponding transition
        U->index = set.size();

        set.insert(U);
        queue.push(std::make_shared<SafraTree>(*U));

        rabin_map.insert({{T->index, i}, U->index});
        dbg(OutputType::General, printf("# f(%u, %u) -> %u\n\n", T->index, i, U->index));
        dbg(OutputType::General, printf("# Unique  %u\n\n", U->index));
      } else {
        // same as a previously computed tree
        // add the corresponding transition
        rabin_map.insert({{T->index, i}, (*itr)->index});
        dbg(OutputType::General, printf("# f(%u, %u) -> %u\n\n", T->index, i, (*itr)->index));
        dbg(OutputType::General, printf("# Index  %u\n\n", (*itr)->index));
      }
    }
  }

  // double check # trees hashed is the same as # trees indexed
  BOOST_ASSERT(num_trees == rabin_map.size()+1);

  if (verbose > static_cast<int>(OutputType::Quiet)) {
    printf("# Trees Generated    %lu\n", num_trees);
    printf("# Trees Hashed       %zd\n", set.size());
    printf("# Empty Trees        %lu\n", num_empty);
    printf("# Rabin Transitions  %zd\n\n", rabin_map.size());
  }

  if (verbose > static_cast<int>(OutputType::General)) {
    printf("# Trees\n");

    for (auto& tree : set) {
      printf("Index: %u\n", tree->index);
      tree->PrintTree();
    }

    // for (auto i = 0UL; i < set.size(); i++) {
    //   for (auto k = 0UL; k < num_alphabet; k++) {
    //     auto itr = rabin_map.find({i, k});
    //     if (itr != rabin_map.end()) {
    //       printf("(%u, %u) -> %u\n", i, k, itr->second);
    //     }
    //   }
    // }

    printf("\n");
  }

  auto rabin = std::make_unique<RabinAutomaton>(num_alphabet, set.size());
  rabin->Init(rabin_map);

  dbg(OutputType::Quiet, printf("\n# PAIRS\n"));

  // Rabin left and right
  // L is set of trees in which state i does not appear
  // R is set of trees in which state i is marked
  boost::dynamic_bitset<> L(set.size());
  boost::dynamic_bitset<> R(set.size());

  for (auto i = 0U; i < 2*num_vertices; i++) {
    L.reset();
    R.reset();

    // Iterate over all trees and fill the bitsets L and R.
    for (auto& tree : set) {
      if (!tree->names[2*i]) {
        L.set(tree->index);
      } else if (tree->names[2*i+1]) {
        R.set(tree->index);
      }
    }

    // L and R must be disjoint.
    BOOST_ASSERT((L & R).none());

    // If R is not empty, then state i appears marked in at least one tree,
    // therefore (L, R) form a valid rabin condition.
    if (R.any()) {
      rabin->pairs.push_back({L, R});

      if (verbose > static_cast<int>(OutputType::Quiet)) {
        std::cout << "i = " << (i/2) << "\n"
                  << rabin->pairs.back().left << "\n"
                  << rabin->pairs.back().right << "\n" << std::endl;
      }
    }
  }

  if (!dotfile.empty()) {
    std::ofstream file(dotfile);
    if (file.fail()) {
      std::cerr << "Failed to open file: " << dotfile << std::endl;
      exit(EXIT_FAILURE);
    }

    boost::write_graphviz(
        file,
        rabin->graph,
        boost::default_writer(),
        label_writer(boost::get(boost::edge_name, rabin->graph)));

    file.close();
  }

  return rabin;
}

boost::dynamic_bitset<> Update(
    boost::dynamic_bitset<>& Q,
    const TransitionMap& map,
    int_type symbol) {
  boost::dynamic_bitset<> L(Q.size());
  //= new boost::dynamic_bitset<>(Q->size());

  // for each state in Q, find all reachable states under the given
  // symbol and add them to the new label set
  ITERATE_BITSET(i, Q) {
    // T.equal_range returns a pair of iterators that give all values that
    // state i maps to under the given symbol
    // if itr == end, then there is no transition
    for (auto [itr, end] = map.equal_range({i, symbol}); itr != end; ++itr) {
      // the second element of the iterator is the new state
      L.set(itr->second);
    }
  }

  return L;
}

std::unique_ptr<BuchiAutomaton> BuchiAutomaton::RabinScott(
    const TransitionMap& map) const {
  TransitionMap dfa_map;

  std::unordered_map<unsigned long, size_type> hash;
  boost::dynamic_bitset<>& Q =
    const_cast<boost::dynamic_bitset<>&>(initial_states);
  boost::dynamic_bitset<> N;

  std::vector<boost::dynamic_bitset<>> state_list;
  std::queue<boost::dynamic_bitset<>> queue;

  size_type num_states = 1;
  size_type num_empty = 0;
  size_type empty_index = 0;
  size_type N_index;

  // push the initial state onto the queue, list, and hash table
  hash.insert({Q.to_ulong(), 0});
  state_list.push_back(Q);
  queue.push(Q);

  while (!queue.empty()) {
    Q = queue.front();

    auto itr = hash.find(Q.to_ulong());

    size_type index = itr->second;

    if (verbose > static_cast<int>(OutputType::General)) {
      printf("Index:  %llu  ", index);
      std::cout << Q << std::endl;
    }

    for (auto i = 0UL; i < num_alphabet; i++) {
      N = Update(Q, map, i);
      num_states++;

      // sink
      if (N.none()) {
        if (num_empty == 0) {
          empty_index = state_list.size();

          hash.insert(std::make_pair(N.to_ulong(), empty_index));
          state_list.push_back(N);

          // The empty tree is an absorbing state with a self-loop under all
          // transition symbols.
          for (auto j = 0UL; j < num_alphabet; j++) {
            dfa_map.insert({{empty_index, j}, empty_index});
          }

          num_empty = num_alphabet;

          if (verbose > static_cast<int>(OutputType::General)) {
            printf("# Unique  %llu  ", empty_index);
            std::cout << N << "\n" << std::endl;
          }
        } else {
          dbg(OutputType::General, printf("\n# Empty  %llu\n\n", empty_index));
        }

        dfa_map.insert({{index, i}, empty_index});
        num_empty++;

        queue.pop();
        continue;
      }

      itr = hash.find(N.to_ulong());

      // unique state
      // - add the bitset to the hash table of states
      // - push the state to the queue of pending states
      // - add the corresponding transition
      if (itr == hash.end()) {
        N_index = state_list.size();

        queue.push(N);
        state_list.push_back(N);

        dfa_map.insert({{index, i}, N_index});

        if (verbose > static_cast<int>(OutputType::General)) {
          printf("# Unique  %llu  ", N_index);
          std::cout << N << "\n\n";
        }

        hash.insert(std::make_pair(N.to_ulong(), N_index));
      } else {
        // same as a previously computed tree
        // add the corresponding transition
        dfa_map.insert({{index, i}, itr->second});
        dbg(OutputType::General, printf("# Index  %llu\n\n", itr->second));
      }
    }

    queue.pop();
  }

  // double check # trees hashed is the same as # trees indexed
  BOOST_ASSERT(state_list.size() == hash.size());

  // container used to store transition function uses an uint32_t to store
  // the key of the tree, so make sure result didn't overflow
  BOOST_ASSERT(state_list.size() <= UINT32_MAX);

  if (verbose > static_cast<int>(OutputType::Quiet)) {
    printf("# STATES GENERATED %llu\n", num_states);
    printf("# STATES HASHED    %zd\n", hash.size());
    printf("# EMPTY STATES     %llu\n", num_empty);
    printf("# TRANSITIONS      %zd\n\n", dfa_map.size());
  }

  dbg(OutputType::General, printf("# STATES\n"));

  // dfa_map.insert({{state_list.size(), 0}, 0});

  auto B = std::make_unique<BuchiAutomaton>(num_alphabet, state_list.size());
  B->Init(dfa_map);

  B->initial_states.set(0);
  B->final_states.reset();

  for (auto i = 0U; i < state_list.size(); i++) {
    boost::dynamic_bitset<> S(state_list[i]);

    if (verbose > static_cast<int>(OutputType::General)) {
      printf("Index: %u  ", i);
      std::cout << S << std::endl;
    }

    if (S.any() && S.intersects(final_states)) {
      B->final_states.set(i);
    }
  }

  dbg(OutputType::General, printf("\n"));

  if (verbose > static_cast<int>(OutputType::Debug)) {
    B->Print();
  }

  return B;
}

void BuchiAutomaton::Minimize() {
  int64_t sink = -1;

  int64_t e0 = -1;
  int64_t e1 = -1;

  std::vector<int> E(num_vertices);

  dbg(OutputType::General, printf("# Minimize()\n\n"));

  char fmt[MAXLINE];
  std::snprintf(fmt, MAXLINE, "%%%dzd  index = %%%dd    ",
      decimal_digits(num_vertices), decimal_digits(num_vertices));

  // Find identical bitsets, assign them to the same class, then proceed with
  // normal minimization algorithm.
  for (auto i = 0U; i < num_vertices; i++) {
    if (!final_states[i]) {
      if (e0 == -1) {
        e0 = i;
        E[i] = i;
        sink = i;
      } else {
        E[i] = e0;
      }
    } else {
      if (e1 == -1) {
        e1 = i;
        E[i] = i;
      } else {
        E[i] = e1;
      }
    }
    dbg(OutputType::General, printf(fmt, i, E[i]));
    dbg(OutputType::General, std::cout << final_states[i] << "\n");
  }
  if (false) {
  std::unordered_map<bool, int> H;
    int i = 0;
    auto H_itr = H.find(final_states[i]);

    if (H_itr != H.end()) {
      E[i] = H_itr->second;

      // if (verbose > OutputType::General) {
      //   printf(fmt, i, E[i]);
      //   std::cout << final_states[i] << "\n";
      // }

      dbg(OutputType::General, printf(fmt, i, E[i]));
      dbg(OutputType::General, std::cout << final_states[i] << "\n");
    } else {
      // assign unique index to the ith pair
      E[i] = i;
      H.insert(std::pair<bool, int>(final_states[i], i));

      if (!final_states[i]) {
        sink = E[i];
      }

      // if (verbose > static_cast<int>(OutputType::General)) {
      //   printf(fmt, i, E[i]);
      //   std::cout << final_states[i] << " (*)\n";
      // }

      dbg(OutputType::General, printf(fmt, i, E[i]));
      dbg(OutputType::General, std::cout << final_states[i] << " (*)\n");
    }
  }

  dbg(OutputType::General, printf("\n"));

  size_type num_rounds = 0;
  size_type num_states = 0;

  snprintf(fmt, MAXLINE, "  %%%dd", std::max(decimal_digits(num_vertices), 2));

  std::vector<int> prev(num_vertices);
  std::vector<int> temp(num_vertices);
  std::vector<int> names(num_vertices);

  auto index = boost::get(boost::vertex_index, graph);
  auto label = boost::get(boost::edge_name, graph);

  bool equivalent = false;
  while (!equivalent) {
    num_rounds++;

    dbg(OutputType::General, printf("%*s:", binary_digits(num_alphabet), "E"));

    // prev is the initial equivalence class at the start of the round
    for (auto i = 0UL; i < num_vertices; i++) {
      temp[i] = -1;
      prev[i] = E[i];

      dbg(OutputType::General, printf(fmt, prev[i]));
    }

    for (auto s = 0UL; s < num_alphabet; s++) {
      std::unordered_map<int_pair, int_type> m(num_vertices);

      dbg(OutputType::General, printf("\n"));
      dbg(OutputType::General, print_binary(s, binary_digits(num_alphabet)));
      dbg(OutputType::General, printf(":"));

      for (auto j = 0UL; j < num_vertices; j++) {
        auto u = boost::vertex(j, graph);

        // find all edges leaving the current vertex and check if any match
        // the current input symbol
        for (auto [itr, end] = boost::out_edges(u, graph); itr != end; ++itr) {
          auto v = boost::target(*itr, graph);

          // if matching edge is found, update class of vertex to target class
          if (s == boost::get(label, *itr)) {
            temp[j] = E[index[v]];
            break;
          }
        }

        dbg(OutputType::General, printf(fmt, temp[j]));
      }

      dbg(OutputType::General, printf("\n%*s:", binary_digits(num_alphabet), "*"));

      size_type current = 0;
      equivalent = true;
      for (auto j = 0UL; j < num_vertices; j++) {
        auto itr = m.find({E[j], temp[j]});

        if (itr == m.end()) {
          m.insert({{E[j], temp[j]}, j});
          E[j] = j;
          names[j] = current++;
        } else {
          E[j] = itr->second;
          names[j] = names[E[j]];
        }

        if (prev[j] != E[j]) {
          equivalent = false;
        }

        temp[j] = -1;
        dbg(OutputType::General, printf(fmt, E[j]));
      }

      num_states = m.size();
      dbg(OutputType::General, printf("\nUnique states = %llu\n", num_states));
    }

    dbg(OutputType::General, printf("\n\n"));
  }

  dbg(OutputType::Quiet, printf("Rounds = %llu\n", num_rounds));

  // once behavioral equivalence is computed, merge equivalent states
  if (num_states == num_vertices) {
    dbg(OutputType::General, printf("All states behaviorally unique.\n\n"));
    return;
  }

  // create a new graph with only the unique states
  // renumber the states so that the unique states are consecutive
  graph_t H(num_states);

  auto label_H = boost::get(boost::edge_name, H);
  for (auto i = 0UL; i < num_vertices; i++) {
    if (E[i] != static_cast<int>(i)) {
      continue;
    }

    auto u = boost::vertex(i, graph);
    for (auto [itr, end] = boost::out_edges(u, graph); itr != end; ++itr) {
      auto v = boost::target(*itr, graph);
      auto j = index[v];

      // u is the source and v is the target
      auto u_H = boost::vertex(names[i], H);
      auto v_H = boost::vertex(names[j], H);

      auto [e, b] = boost::add_edge(u_H, v_H, H);

      // label the transition with the current symbol
      label_H[e] = label[*itr];
    }
  }

  dbg(OutputType::General, printf("States removed = %llu\n\n", num_vertices-num_states));

  graph = H;

  // make the new initial and final state sets
  boost::dynamic_bitset<> I(num_states);
  boost::dynamic_bitset<> F(num_states);

  // I.set(names[initial_states.find_first()]);
  ITERATE_BITSET(i, initial_states) {
    I.set(names[i]);
  }
  // F.set();
  ITERATE_BITSET(i, final_states) {
    F.set(names[i]);
  }

  if (sink != -1) {
    F.reset(names[sink]);
  }

  num_vertices = num_states;
  initial_states = I;
  final_states = F;

  if (verbose > static_cast<int>(OutputType::Quiet)) {
    Print();
  }
}

// Given two Büchi automata, produce the interection automaton.
// Uses a modified pebbling construction similar to intersection of DFAs.
std::unique_ptr<BuchiAutomaton> Intersection(const BuchiAutomaton& A, const BuchiAutomaton& B) {
  auto C = std::make_unique<BuchiAutomaton>(std::min(A.num_alphabet, B.num_alphabet));

  // Speed up graph construction by reserving enough space for vertices to
  // avoid resizing.
  C->Reserve(3*A.num_vertices*B.num_vertices);

  auto label_C = boost::get(boost::edge_name, C->graph);
  auto state_C = boost::get(boost::vertex_name, C->graph);

  auto label_A = boost::get(boost::edge_name, A.graph);
  auto index_A = boost::get(boost::vertex_index, A.graph);

  auto label_B = boost::get(boost::edge_name, B.graph);
  auto index_B = boost::get(boost::vertex_index, B.graph);

  // map: (i, j, component) -> vertex_descriptor
  std::unordered_map<int_triple, graph_t::vertex_descriptor> map;

  std::queue<int_triple> queue;
  std::queue<graph_t::vertex_descriptor> vertex_queue;

  dbg(OutputType::Debug, printf("# N(A) = %llu, N(B) = %llu, S = %llu\n", A.num_vertices, B.num_vertices, C->num_alphabet));
  dbg(OutputType::Debug, printf("# Initial States\n"));

  char fmt_vrtx[MAXLINE];
  snprintf(fmt_vrtx, MAXLINE, "(%%%du, %%%du, ",
           decimal_digits(A.num_vertices), decimal_digits(B.num_vertices));

  ITERATE_BITSET(i, A.initial_states) {
    ITERATE_BITSET(j, B.initial_states) {
      // Add a new vertex to the graph representing the state
      // (i, j, NodeType::Initial).
      auto u_C = boost::add_vertex(C->graph);

      // Store (i, j, NodeType::Initial) -> u in the map.
      auto t = std::make_tuple(i, j, static_cast<int>(NodeType::Initial));
      map.insert({t, u_C});

      // Fill in the graph as a BFS from the initial states.
      queue.push(std::move(t));
      vertex_queue.push(u_C);

      // The state (a, b, NodeType::Initial) is both initial and final.
      state_C[u_C] = NodeType::Initial;

      dbg(OutputType::Debug, printf(fmt_vrtx, i, j));
      dbg(OutputType::Debug, print_state(NodeType::Initial));
      dbg(OutputType::Debug, printf(")\n"));
    }
  }

  dbg(OutputType::Debug, printf("\n"));

  snprintf(fmt_vrtx, MAXLINE, "(%%%du, %%%du, ",
      decimal_digits(A.num_vertices), decimal_digits(B.num_vertices));

  char fmt_edge[MAXLINE];
  snprintf(fmt_edge, MAXLINE, "  -->  (%%%du, %%%du, ",
      decimal_digits(A.num_vertices), decimal_digits(B.num_vertices));

  while (!queue.empty()) {
    auto [i_A, i_B, component] = queue.front();
    auto u_C = vertex_queue.front();
    auto nodetype = static_cast<NodeType>(component);

    queue.pop();
    vertex_queue.pop();

    dbg(OutputType::Debug, printf(fmt_vrtx, i_A, i_B));
    dbg(OutputType::Debug, print_state(static_cast<NodeType>(component)));
    dbg(OutputType::Debug, printf(")\n"));

    auto u_A = boost::vertex(i_A, A.graph);
    auto u_B = boost::vertex(i_B, B.graph);

    if (nodetype == NodeType::Initial) {
      // NodeType::Initial -> NodeType::Final1 on next input
      nodetype = NodeType::Final1;
    } else if (nodetype == NodeType::Final1 && A.final_states[i_A]) {
      // NodeType::Final1 -> FinalType::Final2 when a final state of A is seen
      nodetype = NodeType::Final2;
    } else if (nodetype == NodeType::Final2 && B.final_states[i_B]) {
      // NodeType::Final2 -> NodeType::Initial when a final state of B is seen
      nodetype = NodeType::Initial;
    }

    // Explore out edges of u_A in A and match them with out edges of u_B in B.
    // BGL_FORALL_OUTEDGES(u_A, e_A, A.graph, graph_t) {
    for (auto [itr_A, end_A] = boost::out_edges(u_A, A.graph); itr_A != end_A; ++itr_A) {
      auto e_A = *itr_A;
      auto v_A = boost::target(e_A, A.graph);
      auto i_vA = index_A[v_A];
      auto symbol = label_A[e_A];

      // BGL_FORALL_OUTEDGES(u_B, e_B, B.graph, graph_t) {
      for (auto [itr_B, end_B] = boost::out_edges(u_B, B.graph); itr_B != end_B; ++itr_B) {
        auto e_B = *itr_B;
        if (symbol != label_B[e_B]) {
          continue;
        }

        auto v_B = boost::target(e_B, B.graph);
        auto i_vB = index_B[v_B];

        dbg(OutputType::Debug, printf("  "));
        dbg(OutputType::Debug, print_binary(symbol, binary_digits(C->num_alphabet)));
        dbg(OutputType::Debug, printf(fmt_edge, i_vA, i_vB));
        dbg(OutputType::Debug, print_state(nodetype));
        dbg(OutputType::Debug, printf(")"));

        auto tuple = std::make_tuple(i_vA, i_vB, static_cast<int_type>(nodetype));
        auto itr = map.find(tuple);

        // Add a new vertex to the product machine.
        if (itr == map.end()) {
          auto v_C = boost::add_vertex(C->graph);

          if (nodetype == NodeType::Initial) {
            state_C[v_C] = NodeType::Final;
          } else {
            state_C[v_C] = NodeType::None;
          }

          dbg(OutputType::Debug, printf("  *"));
          itr = map.insert({tuple, v_C}).first;

          queue.push(std::move(tuple));
          vertex_queue.push(v_C);
        }

        // Add a new edge to the target vertex.
        auto [e_C, added] = boost::add_edge(u_C, itr->second, C->graph);
        label_C[e_C] = symbol;

        dbg(OutputType::Debug, printf("\n"));
      }
    }

    dbg(OutputType::Debug, printf("\n"));
  }

  C->Resize();
  dbg(OutputType::Debug, C->Print());

  C->Clean();
  dbg(OutputType::General, C->Print());

  return C;
}

/**
  * Given two Büchi automata, produce the union automaton.
  *
  * Merge the two disjoint graphs into a single graph with two components.
 **/
std::unique_ptr<BuchiAutomaton> DisjointUnion(const BuchiAutomaton& A, const BuchiAutomaton& B) {
  auto num_alphabet = std::max(A.num_alphabet, B.num_alphabet);
  auto num_vertices = A.num_vertices + B.num_vertices;
  auto output = std::make_unique<BuchiAutomaton>(num_alphabet, num_vertices);

  auto label = boost::get(boost::edge_name, output->graph);
  auto state = boost::get(boost::vertex_name, output->graph);

  auto index_A = boost::get(boost::vertex_index, A.graph);
  auto state_A = boost::get(boost::vertex_name, A.graph);
  auto label_A = boost::get(boost::edge_name, A.graph);

  for (auto i = 0UL; i < A.num_vertices; i++) {
    auto u = boost::vertex(i, output->graph);
    auto u_A = boost::vertex(i, A.graph);

    if (A.initial_states[i]) {
      output->initial_states.set(i);
    }

    if (A.final_states[i]) {
      output->final_states.set(i);
    }

    for (auto [itr, end] = boost::out_edges(u_A, A.graph); itr != end; ++itr) {
      auto v_A = boost::target(*itr, A.graph);

      auto index = index_A[v_A];
      auto symbol = label_A[*itr];

      auto v = boost::vertex(index, output->graph);

      auto [e, b] = boost::add_edge(u, v, output->graph);
      label[e] = symbol;
    }

    state[u] = state_A[u_A];
  }

  auto index_B = boost::get(boost::vertex_index, B.graph);
  auto state_B = boost::get(boost::vertex_name, B.graph);
  auto label_B = boost::get(boost::edge_name, B.graph);

  for (auto i = 0UL; i < B.num_vertices; i++) {
    auto u = boost::vertex(i + A.num_vertices, output->graph);
    auto u_B = boost::vertex(i, B.graph);

    if (B.initial_states[i]) {
      output->initial_states.set(i + A.num_vertices);
    }

    if (B.final_states[i]) {
      output->final_states.set(i + A.num_vertices);
    }

    for (auto [itr, end] = boost::out_edges(u_B, B.graph); itr != end; ++itr) {
      auto v_B = boost::target(*itr, B.graph);

      auto index = index_B[v_B] + A.num_vertices;
      auto symbol = label_B[*itr];

      auto v = boost::vertex(index, output->graph);

      auto [e, b] = boost::add_edge(u, v, output->graph);
      label[e] = symbol;
    }

    state[u] = state[u_B];
  }

  output->num_edges = A.num_edges + B.num_edges;

  // may not work
  // output->Clean();
  // output->Resize();

  if (verbose > static_cast<int>(OutputType::General)) {
    output->Print();
  }

  return output;
}

} // namespace omega
