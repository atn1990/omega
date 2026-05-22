#include "Util.h"

#include <boost/dynamic_bitset.hpp>

int verbose = 1;

std::string dotfile = "";
// size_t num_threads = 1;

namespace omega {

std::string print_state(NodeType s) {
  if (s == NodeType::None) {
    return "None";
  } else if (s == NodeType::Initial) {
    return "Initial";
  } else if (s == NodeType::Final) {
    return "Final";
  } else if (s == NodeType::Final1) {
    return "Final1";
  } else if (s == NodeType::Final2) {
    return "Final2";
  } else {
    return "Both";
  }
}

}  // namespace omega
