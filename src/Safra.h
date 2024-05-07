/**
  * <Safra.h>
  *   Class definitions for SafraNode and SafraTree.
 **/

#ifndef OMEGA_SAFRA_H
#define OMEGA_SAFRA_H

#include <cstdint>
#include <iostream>
#include <list>
#include <memory>
#include <ostream>
#include <unordered_map>

#include <boost/dynamic_bitset.hpp>

#include "Automaton.h"

namespace omega {

// using declarations
// single bit, only leaves can be marked
enum class Status : uint8_t {
  Unmarked,
  Marked,
};

class SafraNode {
  public:
    SafraNode(const boost::dynamic_bitset<>&, boost::dynamic_bitset<>&);
    SafraNode(const SafraNode &);

    SafraNode &operator=(const SafraNode&) = delete;
    SafraNode &operator=(SafraNode &&) = delete;

    void Unmark();
    void Update(const TransitionMap& map, int_type symbol);
    void Create(const boost::dynamic_bitset<>&, boost::dynamic_bitset<>&);
    void HorizontalMerge(boost::dynamic_bitset<>& names);
    void KillEmpty();
    void VerticalMerge();

    size_t FillNames(boost::dynamic_bitset<>& names);

    void PrintNode(std::ostream& os, int level) const;

    bool operator==(const SafraNode&) const;
    bool operator!=(const SafraNode&) const;

    friend size_t hash_value(const SafraNode&);

    uint32_t name; // v \in {1, 2, ..., 2n}

    Status status { Status::Unmarked };
    boost::dynamic_bitset<> label; // \emptyset != L(v) \subset Q

    std::list<std::unique_ptr<SafraNode>> children;
};

class SafraTree {
  public:
    SafraTree(int_type num_states, const boost::dynamic_bitset<> &initial_states, const boost::dynamic_bitset<> &final_states);
    SafraTree(const SafraTree& tree);

    SafraTree &operator=(const SafraTree&) = delete;
    SafraTree &operator=(SafraTree&&) = delete;

    void Init(const boost::dynamic_bitset<> &initial_states,
              const boost::dynamic_bitset<> &final_states);

    void Unmark();
    void Update(const TransitionMap& map, int_type symbol);
    void Create(const boost::dynamic_bitset<> &final_states);
    void HorizontalMerge();
    void KillEmpty();
    void VerticalMerge();

    void PrintTree(std::ostream& os = std::cout) const;

    bool operator==(const SafraTree& tree) const;
    friend size_t hash_value(const SafraTree& tree);

    std::unique_ptr<SafraNode> root;

    boost::dynamic_bitset<> names;

    // unique identifier for rabin transitions and rabin pairs
    uint32_t index = 0;

  private:
    uint32_t num_nodes = 0;
};

} // namespace omega

namespace std {
  template <> struct hash<omega::SafraNode> {
    size_t operator()(const omega::SafraNode& node) const;
  };
  template <> struct hash<omega::SafraTree> {
    size_t operator()(const omega::SafraTree& tree) const;
  };
}


#endif // OMEGA_SAFRA_H
