// <Automaton.cpp>
//   Implementation of the Automaton class.

#include "Automaton.h"
#include "Util.h"

#include <cstdio>
#include <iostream>

#include <boost/dynamic_bitset.hpp>

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
  // num_vertices = boost::num_vertices(graph);
  // num_edges = boost::num_edges(graph);
  auto n = decimal_digits(num_vertices);
  auto m = decimal_digits(num_edges);
  auto s = binary_digits(num_alphabet);

  printf("# GRAPH\n");
  printf("  N = %llu\n", num_vertices);
  printf("  M = %llu\n", num_edges);
  printf("  S = %llu\n\n", num_alphabet);

  printf("# ADJACENCY LIST\n");

  char fmt[MAXLINE];
  snprintf(fmt, MAXLINE, "%%%dzd  %%%du  %%%du\n", n, decimal_digits(num_alphabet), n);

  char fmt_vertex[MAXLINE];
  snprintf(fmt_vertex, MAXLINE, "%%%dzd (", n);

  char fmt_v[MAXLINE];
  snprintf(fmt_v, MAXLINE, "  f(%%%dzd, ", n);

  char fmt_edge[MAXLINE];
  snprintf(fmt_edge, MAXLINE, ") = %%%du\n", n);

  // Boost Property Maps
  auto index = boost::get(boost::vertex_index, graph);
  auto label = boost::get(boost::edge_name, graph);
  auto state = boost::get(boost::vertex_name, graph);

  for (auto i = 0U; i < num_vertices; i++) {
    auto u = boost::vertex(i, graph);

    dbg(OutputType::Outfile, printf(fmt_vertex, i));
    print_state(state[u]);
    dbg(OutputType::Outfile, printf(")\n"));

    for (auto [e_itr, e_end] = boost::out_edges(u, graph); e_itr != e_end; ++e_itr) {
      auto v = boost::target(*e_itr, graph);
      auto symbol = label[*e_itr];

      if (verbose == static_cast<int>(OutputType::Outfile)) {
        printf(fmt, i, symbol, index[v]);
      } else {
        printf(fmt_v, index[u]);
        print_binary(symbol, s);
        printf(fmt_edge, index[v]);
        // printf("  ");
        // print_binary(symbol, s);
        // printf(fmt_edge, index[v]);
      }
    }

    printf("\n");
  }

  std::cout << std::endl;
}

} // namespace omega
