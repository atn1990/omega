// <Safra.cpp>
//   Implementation of the SafraNode and SafraTree classes defined in and the
//   core of Safra's determinization algorithm.

#include "Safra.h"
#include "Util.h"

#include <ostream>
#include <memory>
#include <utility>

#include <boost/functional/hash.hpp>

namespace omega {

// SafraNode constructor
// Assigns a node the first unallocated name in the tree.
SafraNode::SafraNode(const boost::dynamic_bitset<>& label, boost::dynamic_bitset<>& names)
  : label(label) {
  for (auto i = 0U; i < names.size(); i += 2) {
    if (!names[i]) {
      name = i/2;
      names.set(i);

      return;
    }
  }

  // If we reach here, no free name was found
  throw std::runtime_error("No free name available for SafraNode");
}


// SafraNode copy constructor
SafraNode::SafraNode(const SafraNode& node) : name(node.name), marked(node.marked), label(node.label) {
  // make a separate copy of each child node
  for (auto& child : node.children) {
    children.emplace_back(std::make_unique<SafraNode>(*child));
  }
}

// Unmarks current node and recursively unmarks its descendants.
void SafraNode::Unmark() {
  marked = false;

  for (auto& child : children) {
    child->Unmark();
  }
}

// Updates the label set according to the transition function under a given
// input symbol.
// The transition function has type T : int * int -> int
void SafraNode::Update(const TransitionMap& map, int_type symbol) {
  boost::dynamic_bitset<> new_label(label.size());

  // for each state in the label set, find all reachable states under the given
  // symbol and add them to the new label set
  ITERATE_BITSET(i, label) {
    // T.equal_range returns a pair of iterators that give all values that
    // state i maps to under the given symbol
    auto range = map.equal_range({i, symbol});

    // if itr.first == itr.second, then there is no transition
    for (auto itr = range.first; itr != range.second; ++itr) {
      // the second element of the iterator is the new state
      new_label.set(itr->second);
    }
  }

  label = new_label;

  // recursively update descendants
  for (auto& child : children) {
    child->Update(map, symbol);
  }
}


// Add a new child if the the label intersects the set of final states.
// Pass the set of names to correctly keep track of which names are in use.
void SafraNode::Create(const boost::dynamic_bitset<>& final_states,
                       boost::dynamic_bitset<>& names) {
  // If a node is marked, then it was created by its parent during this step
  // since all existing nodes in the tree are unmarked at the start, so don't
  // create another child.
  if (marked) {
    return;
  }

  if (label.empty()) {
    return;
  }

  if (label.intersects(final_states)) {
    auto child = std::make_unique<SafraNode>(label & final_states, names);
    child->marked = true;
    names.set(2 * child->name + 1);
    children.emplace_back(std::move(child));
  }

  // Recursively create children for unmarked descendants
  for (auto& child : children) {
    child->Create(final_states, names);
  }
}

// Removes any state that appears to the left of a node.
//
// By assumption, once the parent recursively calls HorizontalMerge() on a
// child, that child's label set is already merged with the left siblings.
//
// The parent ensures that a child's set contains no states outside of its
// own label set.
//
// A kill set is passed to siblings and children that gets updated as the
// algorithm runs. When a parent initially calls HorizontalMerge() on its
// children, the kill set is empty.
//
// As a node finishes cleaning its own state set, it adds its own states to
// the kill set so that its siblings do not keep any of its states, or
// those previously seen.
void SafraNode::HorizontalMerge(boost::dynamic_bitset<>& forbidden) {
  // Remove any states in forbidden since they appear to the left of the
  // tree by assumption.

  // ~forbidden is the set of allowed states
  label &= ~forbidden;

  // At this point, the label set may be empty so let KillEmpty() handle
  // destroying this node and its children.
  if (label.none()) {
    return;
  }

  for (auto& child : children) {
    child->HorizontalMerge(forbidden);
  }

  // add current label to the kill set
  forbidden |= label;
}

// After horizontal merge, some states may be left with empty label sets, so
// delete nodes from the linked list and fix the pointers.
//
// Children are destroyed by the SafraNode destructor so only fix the pointers
// to the siblings.
void SafraNode::KillEmpty() {
  auto itr = children.begin();
  while (itr != children.end()) {
    auto node = itr->get();

    if (node->label.none()) {
      itr = children.erase(itr);
    } else {
      node->KillEmpty();
      ++itr;
    }
  }
}

// Merge child nodes with their parent if the Safra conditions are met.
size_t SafraNode::VerticalMerge(boost::dynamic_bitset<>& names) {
  size_t num_nodes = 1;
  names.set(2*name);
  if (marked) {
    names.set(2*name+1);
  }

  boost::dynamic_bitset<> set(label.size());

  for (auto& child : children) {
    set |= child->label;
  }

  // If the union of the every child label is the parent label, then remove
  // all descendants and mark the parent node.
  if (set == label) {
    marked = true;
    names.set(2*name+1);

    // The destructor will take care of the descendants.
    children.clear();
  }

  for (auto& child : children) {
    num_nodes += child->VerticalMerge(names);
  }

  return num_nodes;
}

void SafraNode::PrintNode(std::ostream& os, int level) const {
  std::string indent(2*level, ' ');
  os << indent << "Name:  " << name << "\n";
  os << indent << "Marked:  " << std::boolalpha << marked << "\n";
  os << indent << "Label:  " << label << '\n';

  if (!children.empty()) {
    os << indent << "Children:\n";
  }

  for (auto& child : children) {
    child->PrintNode(os, level+1);
  }

  os << std::endl;
}

bool SafraNode::operator==(const SafraNode& node) const {
  if (name != node.name) {
    return false;
  }

  if (marked != node.marked) {
    return false;
  }

  if (label != node.label) {
    return false;
  }

  if (children.size() != node.children.size()) {
    return false;
  }

  auto c1 = children.begin();
  auto c2 = node.children.begin();

  while (c1 != children.end() && c2 != node.children.end()) {
    if ((**c1) != (**c2)) {
      return false;
    }

    ++c1;
    ++c2;
  }

  if (c1 != children.end() || c2 != node.children.end()) {
    return false;
  }

  return true;
}

bool SafraNode::operator!=(const SafraNode& other) const {
  return !(*this == other);
}

// Create the names bitset of size 4n, where n is the number of states in the
// Büchi automaton.
//
// For an integer i, names[2i] is set if the name is taken in the tree and
// names[2i+1] is set if the corresponding node is also marked.
SafraTree::SafraTree(
  int_type num_states,
  const boost::dynamic_bitset<>& initial_states,
  const boost::dynamic_bitset<>& final_states)
  : names(4 * num_states) {
  root = std::make_unique<SafraNode>(initial_states, names);
  num_nodes = 1;

  // initialize root according to conditions of Safra's algorithm three cases:
  // - I & F == 0
  // - I & F == I
  // - I & F != 0

  if (initial_states.any() && initial_states.is_subset_of(final_states)) {
    root->marked = true;
    names.set(2 * root->name + 1);
  } else if (initial_states.intersects(final_states)) {
    auto child =
      std::make_unique<SafraNode>(initial_states & final_states, names);
    child->marked = true;
    names.set(2 * child->name + 1);
    root->children.emplace_back(std::move(child));

    num_nodes++;
  }
}

SafraTree::SafraTree(const SafraTree& T) {
  root = std::make_unique<SafraNode>(*T.root);
  names = T.names;
  index = T.index;
  num_nodes = T.num_nodes;
}

SafraTree::SafraTree(SafraTree&& T) : root(std::move(T.root)), names(std::move(T.names)), index(T.index), num_nodes(T.num_nodes) {}

void SafraTree::Unmark() {
  root->Unmark();
}

void SafraTree::Update(const TransitionMap& map, int_type symbol) {
  root->Update(map, symbol);
}

void SafraTree::Create(const boost::dynamic_bitset<>& final_states) {
  root->Create(final_states, names);
}

// Start with the root.
// Initial kill set is empty, which changes accordingly.
void SafraTree::HorizontalMerge() {
  boost::dynamic_bitset<> set(root->label.size());

  root->HorizontalMerge(set);
}

void SafraTree::KillEmpty() {
  if (root->label.none()) {
    root.reset(nullptr);
  } else {
    root->KillEmpty();
  }
}

void SafraTree::VerticalMerge() {
  // once nodes are merged, go through and find all names still present or
  // marked in the tree
  names.reset();
  num_nodes = root->VerticalMerge(names);
}

void SafraTree::PrintTree(std::ostream& os) const {
  os << "Names:  " << names << "\n";
  os << "Nodes:  " << num_nodes << "\n";

  if (root != nullptr) {
    root->PrintNode(os, 1);
  }
}

bool SafraTree::operator==(const SafraTree& tree) const {
  if (num_nodes != tree.num_nodes) {
    return false;
  }

  if (names != tree.names) {
    return false;
  }

  if (root != nullptr && tree.root != nullptr) {
    return *root == *tree.root;
  } else {
    return root == tree.root;
  }
}

bool SafraTree::operator!=(const SafraTree& tree) const {
  return !(*this == tree);
}

size_t SafraNode::hash_value() const {
  size_t seed = 0;

  boost::hash_combine(seed, name);
  boost::hash_combine(seed, static_cast<size_t>(marked));

  std::string str;
  boost::to_string(label, str);
  boost::hash_combine(seed, str);

  for (auto& child : children) {
    boost::hash_combine(seed, child->name);
  }

  return seed;
}

size_t SafraTree::hash_value() const {
  size_t seed = num_nodes;

  std::string str;
  boost::to_string(names, str);
  boost::hash_combine(seed, str);

  if (root != nullptr) {
    boost::hash_combine(seed, root->hash_value());
  }

  return seed;
}

} // namespace omega
