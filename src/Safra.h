// <Safra.h>
//   Class definitions for SafraNode and SafraTree.

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

    size_t hash_value() const;

    // v \in {1, 2, ..., 2n}
    uint32_t name;

    bool marked = false;

    // \emptyset != L(v) \subset Q
    boost::dynamic_bitset<> label;

    std::list<std::unique_ptr<SafraNode>> children;
};

class SafraTree {
  public:
    SafraTree(int_type, const boost::dynamic_bitset<> &, const boost::dynamic_bitset<> &);
    SafraTree(const SafraTree& tree);

    SafraTree &operator=(const SafraTree&) = delete;
    SafraTree &operator=(SafraTree&&) = delete;

    void Unmark();
    void Update(const TransitionMap&, int_type);
    void Create(const boost::dynamic_bitset<> &);
    void HorizontalMerge();
    void KillEmpty();
    void VerticalMerge();

    void PrintTree(std::ostream& os = std::cout) const;

    bool operator==(const SafraTree& tree) const;
    bool operator!=(const SafraTree& tree) const;
    size_t hash_value() const;

    std::unique_ptr<SafraNode> root;

    boost::dynamic_bitset<> names;

    // unique identifier for rabin transitions and rabin pairs
    uint32_t index = 0;
    uint32_t num_nodes = 0;
};

} // namespace omega

#endif // OMEGA_SAFRA_H
