#pragma once

// Util.h
//   Preprocessor constants and macros.

#define BOOST_ENABLE_ASSERT_HANDLER
#define BOOST_ENABLE_ASSERT_DEBUG_HANDLER
#define BOOST_STACKTRACE_LINK

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <ostream>
#include <print>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/assert.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <boost/stacktrace.hpp>

namespace boost {
  inline void assertion_failed_msg(
      char const* expr,
      char const* msg,
      char const* function,
      char const* /*file*/,
      long line) {
    std::cerr << "Expression '" << expr << "' is false in function '"
              << function << "' @" << line << ": "
              << (msg ? msg : "<...>") << ".\nBacktrace:\n"
              << boost::stacktrace::stacktrace() << '\n';
    std::abort();
  }

  inline void assertion_failed(
      char const* expr, char const* function, char const* file, long line) {
    ::boost::assertion_failed_msg(expr, nullptr, function, file, line);
  }
} // namespace boost

namespace std {

template <typename T, typename U>
struct hash<std::pair<T, U>> {
  inline size_t operator()(const std::pair<T, U>& p) const {
    size_t seed = 0;

    boost::hash_combine(seed, p.first);
    boost::hash_combine(seed, p.second);

    return seed;
  }
};

template <typename T, typename U, typename V>
struct hash<std::tuple<T, U, V>> {
  inline size_t operator()(const std::tuple<T, U, V>& p) const {
    size_t seed = 0;

    boost::hash_combine(seed, std::get<0>(p));
    boost::hash_combine(seed, std::get<1>(p));
    boost::hash_combine(seed, std::get<2>(p));

    return seed;
  }
};

} // namespace std

template <>
struct std::formatter<boost::dynamic_bitset<>> : std::formatter<std::string_view> {
    auto format(const boost::dynamic_bitset<>& bs, format_context& ctx) const {
      std::string str;
      boost::to_string(bs, str);
      return std::format_to(ctx.out(), "{}", str);
    }
};

namespace omega {

using int_type = uint64_t;

// The upper and lower (k-1) bits of k-bit numbers define the source and target of de Bruijn map, respectively.
constexpr int_type source(int_type x) { return x >> 1; }
constexpr int_type target(int_type x, unsigned k) { return x & ((1UL << k) - 1); }

// Return the i-th least significant bit of n.
constexpr int_type get_bit(int_type n, unsigned i) { return (n >> i) & 0x1; }

// one-step transition map for ECA: the n-th bit of rule.
constexpr int_type map_bit(int_type rule, unsigned n) { return get_bit(rule, n); }

// combine three bits into one number
constexpr int_type compose_3(int_type x1, int_type x2, int_type x3) {
  return ((x1 & 0x1) << 2) | ((x2 & 0x1) << 1) | (x3 & 0x1);
}

#ifdef DEBUG_BUILD
#  define dbg(N, x) if (verbose > std::to_underlying(N)) { x; }
#else
#  define dbg(N, x)
#endif

inline int decimal_digits(std::size_t n) {
  return n > 1 ? static_cast<int>(std::floor(std::log10(n)) + 1) : 1;
}

inline int binary_digits(std::size_t n) {
  return n > 1 ? static_cast<int>(std::floor(std::log2(n))) + 1 : 1;
}

// Range over the set-bit indices of a dynamic_bitset, for use in a
// range-based for loop: `for (auto i : set_bits(bs)) { ... }`.
class set_bits {
public:
  using size_type = boost::dynamic_bitset<>::size_type;

  explicit set_bits(const boost::dynamic_bitset<>& bs) : bs_(bs) {}

  class iterator {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = size_type;
    using difference_type = std::ptrdiff_t;
    using pointer = const size_type*;
    using reference = size_type;

    iterator(const boost::dynamic_bitset<>& bs, size_type pos)
        : bs_(&bs), pos_(pos) {}

    size_type operator*() const { return pos_; }

    iterator& operator++() {
      pos_ = bs_->find_next(pos_);
      return *this;
    }

    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const iterator& other) const { return pos_ == other.pos_; }
    bool operator!=(const iterator& other) const { return !(*this == other); }

  private:
    const boost::dynamic_bitset<>* bs_;
    size_type pos_;
  };

  iterator begin() const { return iterator(bs_, bs_.find_first()); }
  iterator end() const { return iterator(bs_, boost::dynamic_bitset<>::npos); }

private:
  const boost::dynamic_bitset<>& bs_;
};

#define dbg_var(os, var) \
  (os) << __FILE__ << ":" << __LINE__ << " (" << __func__ << ") " \
       << #var << " = [" << (var) << "]" << std::endl

// Constants
#define MAXLINE 128

#define TRIVIAL -1

enum class OutputType {
  Outfile = 0,
  Quiet,
  General,
  Debug,
};

enum class NodeType {
  None = 0,
  Initial,
  Final,
  Both,
  Final1,
  Final2,
};

// Returns a string representation of the type of state s
std::string_view print_state(NodeType s);

} // namespace omega

// Verbosity level for debugging output.
extern int verbose;

extern std::string dotfile;
// extern size_t num_threads;
