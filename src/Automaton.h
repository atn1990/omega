// Automaton.h
//   Class definition of Automaton.
//   Constants and function-like macros.
//   Data structure definitions.

#ifndef OMEGA_AUTOMATON_H
#define OMEGA_AUTOMATON_H

#include "Util.h"

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

// #include <cstddef>
// #include <iostream>
// #include <iterator>
// #include <memory>
// #include <queue>
// #include <vector>

// Boost C++ headers
// #include <boost/config.hpp>
// #include <boost/bind.hpp>
// #include <boost/range/adaptor/reversed.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/two_bit_color_map.hpp>

namespace omega {

using int_type = uint64_t;
using int_pair = std::pair<int_type, int_type>;
using int_triple = std::tuple<int_type, int_type, int_type>;

using TransitionMap = std::unordered_multimap<int_pair, int_type>;

// Adjacency list graph for representing safra trees as nodes and transitions as labeled edges.
//
// OutEdgeList:     listS
// VertexList:      vecS
// Directed:        bidirectionalS
// VertexProperty:  property<vertex_name_t, uint8_t>
// EdgeProperty:    property<edge_name_t, uint8_t>

using VertexProp = boost::property<boost::vertex_name_t, NodeType>;
using EdgeProp = boost::property<boost::edge_name_t, uint8_t>;
using graph_t =
  boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS,
                        VertexProp, EdgeProp>;

using color_map_t =
  boost::two_bit_color_map<
      boost::property_map<graph_t, boost::vertex_index_t>::const_type>;

using color_t = boost::color_traits<boost::two_bit_color_type>;

class Automaton {
  public:
    using size_type = int_type;

    Automaton(size_type num_alphabet = 0, size_type num_vertices = 0);

    Automaton &operator=(const Automaton &) = delete;
    Automaton &operator=(Automaton &&) = delete;

    void Reserve(size_type num_vertices);
    void Resize();
    void Print() const;

    graph_t graph;
    size_type num_alphabet;
    size_type num_vertices;
    size_type num_edges = 0;
};

} // namespace omega

#endif // OMEGA_AUTOMATON_H
