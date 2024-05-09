#include "Util.h"

#include <boost/dynamic_bitset.hpp>

int verbose = 1;

std::string dotfile = "";
size_t num_threads = 1;

/**
  * Print an integer n in binary with k bits.
 **/
namespace omega {

void print_state(NodeType s, std::ostream& os) {
  if (s == NodeType::None) {
    os << "NONE";
  } else if (s == NodeType::Initial) {
    os << "INITIAL";
  } else if (s == NodeType::Final) {
    os << "FINAL";
  } else {
    os << "BOTH";
  }
}

void print_binary(int_type n, int_type k, std::ostream& os) {
  boost::dynamic_bitset<> set(k, n);
  os << set;
}

}  // namespace omega
