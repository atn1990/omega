#define BOOST_TEST_MODULE OmegaTest

#include <boost/test/included/unit_test.hpp>

#include "BuchiAutomaton.h"
#include "ECA.h"
#include "RabinAutomaton.h"
#include "Util.h"

#include <cstdlib>
#include <utility>

#include <boost/dynamic_bitset.hpp>

using namespace omega;

boost::dynamic_bitset<> make_bitset(std::string s) {
  return boost::dynamic_bitset<>(std::move(s));
}

extern int verbose;
// extern size_t num_threads;

struct GlobalFixture {
  void setup() {
    auto env = std::getenv("VERBOSE");
    if (env != nullptr) {
      char* end;
      verbose = std::strtol(env, &end, 10);
    } else {
      verbose = static_cast<int>(omega::OutputType::General);
    }
    // num_threads = 1;
  }
  void teardown() {}
};

// struct MultiThreaded {
//   void setup() {
//     num_threads = 1;
//   }
//   void teardown() {
//     num_threads = 1;
//   }
// };

BOOST_TEST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_CASE(DeterminizeTest1) {
  BuchiAutomaton B(2, 2);

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
  BOOST_TEST(R->num_edges == 4);

  auto col = { RabinPair("01", "10") };
  BOOST_CHECK_EQUAL_COLLECTIONS(R->pairs.begin(), R->pairs.end(), col.begin(), col.end());
}

BOOST_AUTO_TEST_CASE(DeterminizeTest2) {
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

  B.Print();

  auto R = B.Determinize(m);

  BOOST_TEST(R->num_vertices == 3);
  BOOST_TEST(R->num_edges == 6);

  auto col = { RabinPair("101", "010"), RabinPair("011", "100") };
  BOOST_CHECK_EQUAL_COLLECTIONS(R->pairs.begin(), R->pairs.end(), col.begin(), col.end());
}

BOOST_AUTO_TEST_CASE(DeterminizeTest3) {
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

  B.Print();

  omega::DeterminizeOpts opts;
  opts.SwapUpdateCreate = true;
  auto R = B.Determinize(m, opts);
  R->Minimize();

  BOOST_TEST(R->num_vertices == 4);
  BOOST_TEST(R->num_edges == 12);

  auto col = { RabinPair("0000", "0011") };
  BOOST_CHECK_EQUAL_COLLECTIONS(R->pairs.begin(), R->pairs.end(), col.begin(), col.end());
}

BOOST_AUTO_TEST_CASE(ECATest1, * boost::unit_test::disabled()) {
  Run(110, 1, {});
}

// BOOST_AUTO_TEST_CASE(InjectiveTest) {
//   BOOST_TEST(Injective(110) == false);
// }

// BOOST_AUTO_TEST_CASE(SurjectiveTest) {
//   BOOST_TEST(Surjective(110) == false);
// }

// BOOST_FIXTURE_TEST_CASE(InDegreeTest, MultiThreaded) {
BOOST_AUTO_TEST_CASE(InDegreeTest) {
  BOOST_TEST(InDegree(110, 2) == false);
}

BOOST_AUTO_TEST_CASE(BaseTest) {
  BOOST_TEST(Injective(0) == false);
  BOOST_TEST(Injective(255) == false);
  BOOST_TEST(Injective(110) == false);
  BOOST_TEST(Surjective(0) == false);
  BOOST_TEST(Surjective(255) == false);
  BOOST_TEST(Surjective(110) == false);
  BOOST_TEST(Cycle(0, 1) == true);
  BOOST_TEST(Cycle(255, 1) == true);
  BOOST_TEST(Cycle(110, 1) == true);
  BOOST_TEST(Cycle(0, 2) == false);
  BOOST_TEST(Cycle(255, 2) == false);
  BOOST_TEST(Cycle(110, 2) == true);
  BOOST_TEST(Cycle(110, 3) == true);
  BOOST_TEST(Cycle(110, 5) == true);
  BOOST_TEST(Cycle(110, 7) == true);
  BOOST_TEST(Nilpotent(0, 1) == false);
  BOOST_TEST(Nilpotent(255, 1) == false);
  BOOST_TEST(Nilpotent(204, 1) == true);
  BOOST_TEST(Nilpotent(0, 2) == true);
  BOOST_TEST(Nilpotent(255, 2) == true);
  BOOST_TEST(Nilpotent(110, 2) == false);
  BOOST_TEST(Shift(0, 1, ShiftType::Right) == false);
  BOOST_TEST(Shift(255, 1, ShiftType::Right) == false);
  BOOST_TEST(Shift(16, 1, ShiftType::Right) == true);
  BOOST_TEST(Shift(2, 1, ShiftType::Left) == true);
}
