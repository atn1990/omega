#define BOOST_TEST_MODULE OmegaTest

#include <boost/test/included/unit_test.hpp>

#include "BuchiAutomaton.h"
#include "ECA.h"
#include "RabinAutomaton.h"
#include "Util.h"

#include <cstdlib>
#include <format>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include <boost/dynamic_bitset.hpp>

using namespace omega;

boost::dynamic_bitset<> make_bitset(std::string s) {
  return boost::dynamic_bitset<>(std::move(s));
}

extern int verbose;

struct GlobalFixture {
  void setup() {
    auto env = std::getenv("VERBOSE");
    if (env != nullptr) {
      char* end;
      verbose = std::strtol(env, &end, 10);
    } else {
      verbose = std::to_underlying(omega::OutputType::General);
    }
  }
  void teardown() {}
};

BOOST_TEST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_SUITE(DeterminizeTests)

BOOST_AUTO_TEST_CASE(DeterminizeTest1) {
  BuchiAutomaton B(2, 2);

  // Accepts strings with a finite number of 1s followed by infinitely many 0s
  TransitionMap m = {
    {{0, 0}, 0},
    {{0, 1}, 0},
    {{0, 0}, 1},
    {{1, 0}, 1},
  };

  B.initial_states = make_bitset("01");
  B.final_states = make_bitset("10");
  B.Init(m);

  B.Print();

  auto R = B.Determinize(m);

  BOOST_CHECK_EQUAL(R->num_vertices, 2);
  BOOST_CHECK_EQUAL(R->num_edges, 4);

  auto col = { RabinPair("01", "10") };
  BOOST_CHECK_EQUAL_COLLECTIONS(R->pairs.begin(), R->pairs.end(), col.begin(), col.end());
}

BOOST_AUTO_TEST_CASE(DeterminizeTest2) {
  BuchiAutomaton B(2, 2);

  // Accepts strings with a finite number of 1s (at least 1) followed by infinitely many 0s
  TransitionMap m = {
    {{0, 0}, 0},
    {{0, 1}, 0},
    {{0, 1}, 1},
    {{1, 0}, 1},
  };

  B.initial_states = make_bitset("01");
  B.final_states = make_bitset("10");
  B.Init(m);

  B.Print();

  auto R = B.Determinize(m);

  BOOST_CHECK_EQUAL(R->num_vertices, 3);
  BOOST_CHECK_EQUAL(R->num_edges, 6);

  auto col = { RabinPair("101", "010"), RabinPair("011", "100") };
  BOOST_CHECK_EQUAL_COLLECTIONS(R->pairs.begin(), R->pairs.end(), col.begin(), col.end());
}

BOOST_AUTO_TEST_CASE(DeterminizeTest3) {
  BuchiAutomaton B(3, 2);

  // Accepts strings where every occurrence of the symbol 2 is followed, later in the word, by at least one 0
  TransitionMap m = {
    {{0, 0}, 0},
    {{0, 1}, 0},
    {{0, 1}, 1},
    {{0, 2}, 1},
    {{1, 0}, 0},
    {{1, 1}, 1},
    {{1, 2}, 1},
  };

  B.initial_states = make_bitset("01");
  B.final_states = make_bitset("01");
  B.Init(m);

  B.Print();

  omega::DeterminizeOpts opts;
  opts.SwapUpdateCreate = true;
  auto R = B.Determinize(m, opts);
  R->Minimize();

  BOOST_CHECK_EQUAL(R->num_vertices, 4);
  BOOST_CHECK_EQUAL(R->num_edges, 12);

  auto col = { RabinPair("0000", "0011") };
  BOOST_CHECK_EQUAL_COLLECTIONS(R->pairs.begin(), R->pairs.end(), col.begin(), col.end());
}

BOOST_AUTO_TEST_SUITE_END()

// --- New test coverage ---------------------------------------------------
//
// Tests below are grouped by topic. Whenever an assertion is a numeric
// quantity that would normally be hard to hand-derive (num_vertices,
// num_edges, etc. for new automata) it is justified by an invariant that
// must hold for any correct implementation -- never by a hand-computed
// determinization result. The few places where we compare against a
// concrete "characterization" value follow the same style as the existing
// DeterminizeTest1/2/3 cases.

BOOST_AUTO_TEST_SUITE(DeterminizeInvariants)

// Determinization is a deterministic function of its input: running it
// twice on the same input must produce identical Rabin automata
// (same vertex/edge counts and the same set of accepting pairs).
BOOST_AUTO_TEST_CASE(DeterminizeIsReproducible) {
  BuchiAutomaton B(2, 2);
  TransitionMap m = {
    {{0, 0}, 0},
    {{0, 1}, 0},
    {{0, 1}, 1},
    {{1, 0}, 1},
  };
  B.initial_states = make_bitset("01");
  B.final_states = make_bitset("10");
  B.Init(m);

  auto R1 = B.Determinize(m);

  BuchiAutomaton B2(2, 2);
  B2.initial_states = make_bitset("01");
  B2.final_states = make_bitset("10");
  B2.Init(m);
  auto R2 = B2.Determinize(m);

  BOOST_TEST(R1->num_vertices == R2->num_vertices);
  BOOST_TEST(R1->num_edges == R2->num_edges);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      R1->pairs.begin(), R1->pairs.end(),
      R2->pairs.begin(), R2->pairs.end());
}

// Every RabinPair produced by Determinize must have left/right bitsets
// whose width matches the determinized state space (num_vertices), and
// all pairs in the same automaton must share that width. Also, the
// resulting automaton must be non-trivial (>0 vertices, >0 edges).
BOOST_AUTO_TEST_CASE(DeterminizePairsHaveUniformWidth) {
  BuchiAutomaton B(2, 2);
  TransitionMap m = {
    {{0, 0}, 0},
    {{0, 1}, 0},
    {{0, 1}, 1},
    {{1, 0}, 1},
  };
  B.initial_states = make_bitset("01");
  B.final_states = make_bitset("10");
  B.Init(m);

  auto R = B.Determinize(m);

  BOOST_TEST(R->num_vertices > 0u);
  BOOST_TEST(R->num_edges > 0u);
  BOOST_TEST(R->pairs.size() > 0u);

  for (const auto& p : R->pairs) {
    BOOST_TEST(p.left.size() == R->num_vertices);
    BOOST_TEST(p.right.size() == R->num_vertices);
  }
}

// Minimize() should be idempotent: minimizing an already-minimized
// automaton must not change its size or accepting pairs.
BOOST_AUTO_TEST_CASE(MinimizeIsIdempotent) {
  BuchiAutomaton B(3, 2);
  TransitionMap m = {
    {{0, 0}, 0},
    {{0, 1}, 0},
    {{0, 1}, 1},
    {{0, 2}, 1},
    {{1, 0}, 0},
    {{1, 1}, 1},
    {{1, 2}, 1},
  };
  B.initial_states = make_bitset("01");
  B.final_states = make_bitset("01");
  B.Init(m);

  omega::DeterminizeOpts opts;
  opts.SwapUpdateCreate = true;
  auto R = B.Determinize(m, opts);
  R->Minimize();

  const auto v1 = R->num_vertices;
  const auto e1 = R->num_edges;
  const auto pairs1 = R->pairs;

  R->Minimize();

  BOOST_TEST(R->num_vertices == v1);
  BOOST_TEST(R->num_edges == e1);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      R->pairs.begin(), R->pairs.end(),
      pairs1.begin(), pairs1.end());
}

BOOST_AUTO_TEST_SUITE_END()

// -------------------------------------------------------------------------
// RabinPair
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(RabinPairTests)

// The (size_type) ctor must produce two zero-initialized bitsets of the
// requested width.
BOOST_AUTO_TEST_CASE(NumBitsCtorIsZeroInitialized) {
  RabinPair p(5);
  BOOST_TEST(p.left.size() == 5u);
  BOOST_TEST(p.right.size() == 5u);
  BOOST_TEST(p.left.none());
  BOOST_TEST(p.right.none());
}

// The (string,string) and (bitset,bitset) ctors must agree, since
// dynamic_bitset(std::string) is the documented way to build one from a
// 0/1 string.
BOOST_AUTO_TEST_CASE(StringAndBitsetCtorsAgree) {
  RabinPair a("101", "010");
  RabinPair b(make_bitset("101"), make_bitset("010"));
  BOOST_TEST(a == b);
  BOOST_TEST(!(a != b));
}

// Equality is reflexive and symmetric; inequality is its negation.
BOOST_AUTO_TEST_CASE(EqualityIsReflexiveAndSymmetric) {
  RabinPair a("101", "010");
  RabinPair b("101", "010");
  RabinPair c("100", "010");

  BOOST_TEST(a == a);
  BOOST_TEST(a == b);
  BOOST_TEST(b == a);
  BOOST_TEST(a != c);
  BOOST_TEST(c != a);
  BOOST_TEST(!(a == c));
}

// Reset() clears both sides while preserving width.
BOOST_AUTO_TEST_CASE(ResetClearsBitsKeepsWidth) {
  RabinPair p("101", "110");
  const auto w_left = p.left.size();
  const auto w_right = p.right.size();

  p.Reset();

  BOOST_TEST(p.left.size() == w_left);
  BOOST_TEST(p.right.size() == w_right);
  BOOST_TEST(p.left.none());
  BOOST_TEST(p.right.none());
}

// Resize() resets bits and changes both sides to the requested width.
BOOST_AUTO_TEST_CASE(ResizeChangesWidthAndResets) {
  RabinPair p("1010", "0101");
  p.Resize(7);

  BOOST_TEST(p.left.size() == 7u);
  BOOST_TEST(p.right.size() == 7u);
  BOOST_TEST(p.left.none());
  BOOST_TEST(p.right.none());
}

// operator<< must round-trip through std::ostringstream and produce a
// non-empty string for non-trivial inputs.
BOOST_AUTO_TEST_CASE(StreamOperatorProducesOutput) {
  RabinPair p("101", "010");
  std::ostringstream os;
  os << p;
  const auto s = os.str();
  BOOST_TEST(!s.empty());
}

BOOST_AUTO_TEST_SUITE_END()

// -------------------------------------------------------------------------
// std::hash specializations and std::formatter for dynamic_bitset
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(UtilTests)

// Equal keys must hash equal (the basic hash invariant). Combined with
// operator==, that's enough for use as an unordered_map key.
BOOST_AUTO_TEST_CASE(PairHashIsConsistent) {
  std::hash<std::pair<int, int>> h;
  std::pair<int, int> a{3, 7};
  std::pair<int, int> b{3, 7};

  BOOST_TEST(h(a) == h(b));

  std::unordered_map<std::pair<int, int>, std::string> m;
  m[{1, 2}] = "x";
  m[{3, 4}] = "y";
  BOOST_TEST(m.at({1, 2}) == "x");
  BOOST_TEST(m.at({3, 4}) == "y");
}

BOOST_AUTO_TEST_CASE(TupleHashIsConsistent) {
  std::hash<std::tuple<int, int, int>> h;
  std::tuple<int, int, int> a{1, 2, 3};
  std::tuple<int, int, int> b{1, 2, 3};

  BOOST_TEST(h(a) == h(b));

  std::unordered_map<std::tuple<int, int, int>, int> m;
  m[{1, 2, 3}] = 10;
  m[{4, 5, 6}] = 20;
  BOOST_TEST((m.at({1, 2, 3}) == 10));
  BOOST_TEST((m.at({4, 5, 6}) == 20));
}

// std::format("{}", bs) must agree with boost::to_string(bs) since that
// is exactly what the formatter delegates to.
BOOST_AUTO_TEST_CASE(BitsetFormatterMatchesBoostToString) {
  auto bs = make_bitset("10110");
  std::string expected;
  boost::to_string(bs, expected);
  BOOST_TEST(std::format("{}", bs) == expected);

  auto empty = boost::dynamic_bitset<>{};
  std::string expected_empty;
  boost::to_string(empty, expected_empty);
  BOOST_TEST(std::format("{}", empty) == expected_empty);
}

// print_state must produce the exact strings the implementation defines
// for every NodeType enumerator.
BOOST_AUTO_TEST_CASE(PrintStateCoversAllNodeTypes) {
  BOOST_TEST(print_state(NodeType::None) == "None");
  BOOST_TEST(print_state(NodeType::Initial) == "Initial");
  BOOST_TEST(print_state(NodeType::Final) == "Final");
  BOOST_TEST(print_state(NodeType::Both) == "Both");
  BOOST_TEST(print_state(NodeType::Final1) == "Final1");
  BOOST_TEST(print_state(NodeType::Final2) == "Final2");
}

BOOST_AUTO_TEST_SUITE_END()

// -------------------------------------------------------------------------
// ECA predicates
// -------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(ECA_Predicates)

// Surjectivity / injectivity for the trivial constant rules: rule 0 maps
// everything to all-zeros, rule 255 maps everything to all-ones, so
// neither is injective nor surjective
// Rule 204 is the identity, so it is both (bijective)
// Rule 51 is the negation of the identity, so it is also both
BOOST_AUTO_TEST_CASE(InjectiveSurjectiveTrivialRules) {
  for (auto rule : {0, 255}) {
    BOOST_TEST(Injective(rule) == false);
    BOOST_TEST(Surjective(rule) == false);
  }

  for (auto rule : {51, 204}) {
    BOOST_TEST(Injective(rule) == true);
    BOOST_TEST(Surjective(rule) == true);
  }
}

// The all-zero configuration is a 1-cycle of every rule whose lookup table
// sends 000 -> 0 (i.e. rule bit 0 is 0). Cycle(rule, 1) is equivalent to
// FixedPoint(rule) (the implementation imposes no inequality constraint when
// k == 1), so we restate that link here
// Rule 51 sets bit 0 = 1, so 000 -> 1 and the all-zero configuration is *not*
// fixed; in fact rule 51 has no fixed point at all, so Cycle(51, 1) is false
// Rule 0 sends every cell to 0, so the all-zero configuration is fixed
// Rule 255 sends every cell to 1, so the all-one configuration is fixed
// Rule 204 is the identity rule (new cell == middle cell), so every
// configuration is fixed
// Rule 51 has lookup table 0b00110011: bit i of the rule selects new cell
// value for neighborhood i (= 4*l + 2*m + r)
// The output equals (NOT m), independent of l and r
// A configuration x is a fixed point iff x_i == NOT x_i for every i, which is
// impossible -- therefore no fixed point exists
BOOST_AUTO_TEST_CASE(CycleOneMatchesFixedPoint) {
  std::map<int_type, bool> expected = {
    {0, true},
    {51, false},
    {204, true},
    {255, true},
  };

  for (auto [rule, result] : expected) {
    BOOST_TEST(Cycle(rule, 1) == result);
    BOOST_CHECK_EQUAL(Cycle(rule, 1), FixedPoint(rule));
  }
}

BOOST_AUTO_TEST_CASE(ECATest1, * boost::unit_test::disabled()) {
  Run(110, 1, {});
}

BOOST_AUTO_TEST_CASE(BaseTest) {
  BOOST_TEST(Injective(110) == false);
  BOOST_TEST(Surjective(110) == false);
  BOOST_TEST(Cycle(110, 1) == true);
  BOOST_TEST(Cycle(0, 2) == false);
  BOOST_TEST(Cycle(255, 2) == false);
  BOOST_TEST(Cycle(110, 2) == true);
  BOOST_TEST(Cycle(110, 3) == true);
  BOOST_TEST(Cycle(110, 5) == true);
  BOOST_TEST(Cycle(110, 7) == true);
  BOOST_TEST(Nilpotent(0, 1) == true);
  BOOST_TEST(Nilpotent(255, 1) == true);
  BOOST_TEST(Nilpotent(204, 1) == true);
  BOOST_TEST(Nilpotent(110, 2) == false);
  BOOST_TEST(Shift(0, 1, ShiftType::Right) == false);
  BOOST_TEST(Shift(255, 1, ShiftType::Right) == false);
  BOOST_TEST(Shift(16, 1, ShiftType::Right) == true);
  BOOST_TEST(Shift(2, 1, ShiftType::Left) == true);
}


BOOST_AUTO_TEST_SUITE_END()
