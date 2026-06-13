#include "Util.h"

#include <boost/dynamic_bitset.hpp>

int verbose = 1;

std::string dotfile = "";
// size_t num_threads = 1;

namespace omega {

std::string_view print_state(NodeType s) {
  switch (s) {
    case NodeType::None:
      return "None";
    case NodeType::Initial:
      return "Initial";
    case NodeType::Final:
      return "Final";
    case NodeType::Final1:
      return "Final1";
    case NodeType::Final2:
      return "Final2";
    case NodeType::Both:
      return "Both";
  }
}

}  // namespace omega
