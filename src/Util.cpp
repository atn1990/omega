#include "Util.h"

#include <boost/dynamic_bitset.hpp>

int verbose = 1;

std::string dotfile = "";
size_t num_threads = 1;

namespace omega {

void print_state(NodeType s, std::ostream& os) {
  if (s == NodeType::None) {
    os << "None";
  } else if (s == NodeType::Initial) {
    os << "Initial";
  } else if (s == NodeType::Final) {
    os << "Final";
  } else {
    os << "Both";
  }
}

void print_binary(int_type n, int_type k, std::ostream& os) {
  boost::dynamic_bitset<> set(k, n);
  os << set;
}

}  // namespace omega
