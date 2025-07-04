// Util.h
//   Preprocessor constants and macros.

#ifndef OMEGA_UTIL_H
#define OMEGA_UTIL_H

#define BOOST_ENABLE_ASSERT_HANDLER
#define BOOST_ENABLE_ASSERT_DEBUG_HANDLER
#define BOOST_STACKTRACE_LINK

#include <cstdint>
#include <iostream>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>

#include <boost/assert.hpp>
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

namespace omega {

using int_type = uint64_t;

// The upper and lower (k-1) bits of k-bit numbers define the source and target of de Bruijn map, respectively.
#define SOURCE(x) ((x) >> 1)
#define TARGET(x, k) ((x) & ((1 << (2 * (k))) - 1))

// one-step transition map for ECA
#define MAP(rule, n) (static_cast<int_type>(((rule) >> (n)) & 0x1))

// Return the i-th least significant bit of n
#define GET_BIT(n, i) (static_cast<int_type>(((n) >> (i)) & 0x1))

// combine three bits into one number
// #define Q3(x1, x2, x3) (((x1) << 2) | ((x2) << 1) | (x3))
#define COMPOSE_3(x1, x2, x3) \
  ((((x1) & 0x1) << 2) | (((x2) & 0x1) << 1) | ((x3) & 0x1))

// simple conditional logging
// #ifdef DEBUG_BUILD
// #  define dbg_printf(N, ...) if (verbose > N) { printf(__VA_ARGS__); }
// #else
// #  define dbg_printf(N, ...)
// #endif

// #ifdef DEBUG_BUILD
// #  define dbg_cout(N, x) if (verbose > N) { std::cout << x; }
// #else
// #  define dbg_cout(N, x)
// #endif

#ifdef DEBUG_BUILD
#  define dbg(N, x) if (verbose > static_cast<int>(N)) { x; }
#else
#  define dbg(N, x)
#endif

#define decimal_digits(n) \
  ((n) > 1 ? static_cast<int>(std::ceil(std::log10(n))) : 1)

#define binary_digits(n) \
  ((n) > 1 ? static_cast<int>(std::ceil(std::log2(n))) : 1)

#define ITERATE_BITSET(var, bitset) \
  for (auto (var) = (bitset).find_first(); (var) != (bitset).npos; (var) = (bitset).find_next(var))

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

// Print the type of state s to the output stream os.
void print_state(NodeType s, std::ostream& os = std::cout);
// void print_state(int_type s, std::ostream& os = std::cout);

// Print an integer n in binary with k bits.
void print_binary(int_type n, int_type k, std::ostream& os = std::cout);

} // namespace omega

// Verbosity level for debugging output.
extern int verbose;

extern std::string dotfile;
extern size_t num_threads;

#endif // OMEGA_UTIL_H
