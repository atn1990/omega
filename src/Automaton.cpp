// <Automaton.cpp>
//   Implementation of the Automaton class.

#include "Automaton.h"
#include "Util.h"

#include <cstdio>
#include <print>

namespace omega {

Automaton::Automaton(size_type num_alphabet, size_type num_vertices)
  : graph(num_vertices), num_alphabet(num_alphabet), num_vertices(num_vertices) {}

void Automaton::Reserve(size_type num_vertices) {
  graph.m_vertices.reserve(num_vertices);
}

void Automaton::Resize() {
  num_vertices = boost::num_vertices(graph);
  num_edges = boost::num_edges(graph);
}

// Pretty printing of automata.
void Automaton::Print(void) const {
  auto vertex_width = decimal_digits(num_vertices-1);
  auto binary_width = binary_digits(num_alphabet-1);
  auto decimal_width = decimal_digits(num_alphabet-1);

  std::print("# GRAPH\n");
  std::print("  N = {}\n", num_vertices);
  std::print("  M = {}\n", num_edges);
  std::print("  S = {}\n\n", num_alphabet);

  std::print("# ADJACENCY LIST\n");

  // Boost Property Maps
  auto index = boost::get(boost::vertex_index, graph);
  auto label = boost::get(boost::edge_name, graph);
  auto state = boost::get(boost::vertex_name, graph);

  for (auto i = 0U; i < num_vertices; i++) {
    auto u = boost::vertex(i, graph);

    dbg(OutputType::Outfile, std::print("{:{}} ({})\n", i, vertex_width, print_state(state[u])));

    for (auto [e_itr, e_end] = boost::out_edges(u, graph); e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, graph);
      auto symbol = label[*e_itr];

      if (verbose == std::to_underlying(OutputType::Outfile)) {
        std::print("{:{}}  {:{}}  {:{}}", i, vertex_width, symbol, decimal_width, index[v], vertex_width);
      } else {
        std::print("  f({:{}}, {:0{}b}) = {:{}}\n", index[u], vertex_width, symbol, binary_width, index[v], vertex_width);
      }
    }

    std::print("\n");
  }

  std::print("\n");
}

} // namespace omega
