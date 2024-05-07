#ifndef OMEGA_RABIN_AUTOMATON_H
#define OMEGA_RABIN_AUTOMATON_H

#include <cstdint>
#include <list>
#include <string>
#include <unordered_map>
#include <utility>

#include <boost/dynamic_bitset.hpp>

#include "Automaton.h"

namespace omega {

// Deterministic transition system
using RabinTransitionMap = std::unordered_map<int_pair, int_type>;

using RabinPair =
struct RabinPair {
  boost::dynamic_bitset<> left;
  boost::dynamic_bitset<> right;

  using size_type = boost::dynamic_bitset<>::size_type;

  RabinPair() = default;
  explicit RabinPair(size_type num_bits) : left(num_bits), right(num_bits) {}
  RabinPair(boost::dynamic_bitset<> l, boost::dynamic_bitset<> r)
    : left(std::move(l)), right(std::move(r)) {}
  RabinPair(std::string l, std::string r)
    : left(std::move(l)), right(std::move(r)) {}

  void Reset() {
    left.reset();
    right.reset();
  }

  void Resize(size_type num_bits) {
    Reset();
    left.resize(num_bits);
    right.reset(num_bits);
  }

  bool operator==(const RabinPair& p) const {
    return left == p.left && right == p.right;
  }
  bool operator!=(const RabinPair& p) const {
    return !(*this == p);
  }
};

std::ostream& operator<<(std::ostream& os, const RabinPair& p);

class RabinAutomaton : public Automaton {
  public:
    std::list<RabinPair> pairs;

    RabinAutomaton(size_type alphabet = 0, size_type vertices = 0);

    void Init(const RabinTransitionMap& map);

    bool Universal();
    void Clean();
    void Minimize();

    void TestUV(const std::string&, const std::string&);
    void Print() const;
};

} // namespace omega

#endif // OMEGA_RABIN_AUTOMATON_H
