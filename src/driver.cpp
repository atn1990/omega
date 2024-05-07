/**
  * <driver.cpp>
  *   Facilities for reading automata from files and processing command-line
  *   options.
 **/

#include "BuchiAutomaton.h"
#include "RabinAutomaton.h"
#include "ECA.h"
#include "Util.h"

#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
// #include <string_view>
// #include <sstream>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>

// #include <gperftools/profiler.h>

extern int verbose; // sets verbosity level for debugging output
extern size_t num_threads; // number of threads to execute concurrently
extern std::string dotfile;

using omega::BuchiAutomaton;
using omega::RabinAutomaton;
using omega::RabinTransitionMap;
using omega::TransitionMap;

/**
  * Read a "u:v" pair from infile (may be stdin) and parse the input.
 **/
bool ReadUV(std::istream& file,
            const uint32_t alphabet,
            std::string& U,
            std::string& V) {
  if (!std::getline(file, U, ':') || !std::getline(file, V)) {
    return false;
  }

  auto pred =
    [alphabet](unsigned char c) {
      return std::isdigit(c) && (c - '0') < alphabet;
    };

  if (!std::all_of(U.begin(), U.end(), pred) ||
      !std::all_of(V.begin(), V.end(), pred)) {
    V.clear();
  }

  return true;
}


/**
  * Output a classifcation of the ECA from an input trace file.
 **/
void classify(const std::string& filename) {
  // maps a trace to the rules having that trace
  std::unordered_map<std::string, std::list<uint8_t>> map;

  // sort traces according to the lowest numbered rule in its class
  /*
  auto cmp =
    [&](const std::string &tr1, const std::string &tr2) {
      return map.find(tr1)->second.front() < map.find(tr2)->second.front();
    };

  decltype(cmp)
  */

  // maps a size of a class to a set of traces of that class size
  std::map<uint8_t, std::set<std::string>> class_map;

  std::ifstream file(filename);
  if (file.fail()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string str;
  while (std::getline(file, str)) {
    if (str.front() == '#') {
      std::cout << str << std::endl;
      continue;
    }

    //std::istringstream stream(str);
    const char *ptr = str.c_str();
    char *endptr = nullptr;

    uint8_t rule = std::strtol(ptr, &endptr, 10);
    printf("%3u: ", rule);
    // auto start = str.begin();
    // std::advance(start, endptr - ptr);
    // str.erase(str.begin(), start);

    str = endptr;

    // std::isspace is overloaded
    auto pred = [](unsigned char c) { return std::isspace(c); };

    // remove whitespace from str
    str.erase(std::remove_if(str.begin(), str.end(), pred), str.end());
    std::cout << str << std::endl;

    auto itr = map.find(str);

    // if the trace is not found, create a new empty list
    if (itr == map.end()) {
      auto [itr, added] = map.emplace(str, std::list<uint8_t>{});
      BOOST_ASSERT(added);
      // auto [itr, added] = map.emplace(str, std::list<uint8_t>{});
    }

    (void) itr->second.push_back(rule);
  }

  printf("\n\n");
  for (auto trace : map) {
    auto itr = class_map.find(trace.second.size());

    // if the size is not found, create a new empty set
    if (itr == class_map.end()) {
      auto [itr, added] = class_map.emplace(trace.second.size(), std::set<std::string>{});
      BOOST_ASSERT(added);
    }

    (void) itr->second.insert(trace.first);
  }

  // the map stores keys in ascending order, so traverse the class_map in
  // reverse to output the largest classes first
  for (auto r_itr = class_map.rbegin(); r_itr != class_map.rend(); ++r_itr) {
    printf("Class: %u\nSize: %lu\n", r_itr->first, r_itr->second.size());

    for (auto c : r_itr->second) {
      for (auto n : map.find(c)->second) {
        printf("%3d  ", n);
      }

      printf("\n\n");
    }

    printf("\n");
  }

  file.close();
}

/**
  * Convert description of a Büchi automaton from a file.
 **/
std::unique_ptr<BuchiAutomaton> ReadBuchi(const std::string& filename, TransitionMap& map) {
  std::ifstream file(filename);
  if (file.fail()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  uint32_t alphabet = 0;
  uint32_t states = 0;
  uint32_t transitions = 0;
  uint32_t count = 0;

  boost::dynamic_bitset<> I;
  boost::dynamic_bitset<> F;

  std::string str;
  while (std::getline(file, str)) {
    if (str.front() != '#') {
      continue;
    }

    if (str == "# STATES") {
      file >> states;
      std::getline(file, str);
    } else if (str == "# ALPHABET") {
      file >> alphabet;
      std::getline(file, str);

      // optimization: use uint8_t to represent transition symbols on edges
      assert(alphabet <= UINT8_MAX);
    } else if (str == "# TRANSITIONS") {
      file >> transitions;
      std::getline(file, str);
    } else if (str == "# BEGIN TRANSITIONS") {
      uint32_t u, v, s;

      dbg(GENERAL, printf("# TRANSITIONS\n"));

      //std::getline(file, str);
      while (std::getline(file, str) && str != "# END TRANSITIONS") {
        // ignore empty lines
        if (str.empty()) {
          continue;
        }

        if (sscanf(str.c_str(), "%u %u %u", &u, &s, &v) != 3) {
          std::cerr << "Malformed transition: " << str << std::endl;
          exit(EXIT_FAILURE);
        }

        map.insert({{u, s}, v});
        if (verbose > GENERAL) {
          std::cout << "(" << u << ", " << s << ") -> " << v << "\n";
        }

        count++;
      }

      dbg(GENERAL, printf("\n"));
    } else if (str == "# INITIAL") {
      file >> I;
      std::getline(file, str);

      if (verbose > GENERAL) {
        std::cout << "# INITIAL\n" << I << "\n\n";
      }
    } else if (str == "# FINAL") {
      file >> F;
      std::getline(file, str);

      if (verbose > GENERAL) {
        std::cout << "# FINAL\n" << F << "\n\n";
      }
    }
  }

  std::flush(std::cout);

  file.close();

  // Remove useless final states. If they're all final, don't bother.
  if (F.count() != F.size()) {
    dbg(GENERAL, printf("# Removing useless final states\n\n"));
    assert(F.count() != 0);
  }

  assert(states != 0);
  assert(alphabet != 0);
  assert(transitions != 0);
  assert(count == transitions);
  assert(I.size() == states);
  assert(F.size() == states);
  assert(!I.empty());
  assert(!F.empty());

  if (verbose > GENERAL) {
    printf("# STATES      %u\n", states);
    printf("# ALPHABET    %u\n", alphabet);
    printf("# TRANSITIONS %u\n\n", transitions);
  }

  std::unique_ptr<BuchiAutomaton> B = std::make_unique<BuchiAutomaton>(alphabet, states);
  B->Init(map);
  B->initial_states = I;
  B->final_states = F;

  if (verbose > GENERAL) {
    B->Print();
  }

  return B;
}

/**
  * Read two Büchi automata from a file and construct the product automaton.
 **/
void Intersection(const std::vector<std::string> &input) {
  TransitionMap map_A, map_B;

  auto A = ReadBuchi(input.front(), map_A);
  auto B = ReadBuchi(input.back(), map_B);
  auto output = Intersection(*A, *B);

  if (verbose != QUIET) {
    output->Print();
  }
}

/**
  * Read two Büchi automata from a file and construct the union automaton.
 **/
void DisjointUnion(const std::vector<std::string> &input) {
  TransitionMap map_A, map_B;

  auto A = ReadBuchi(input.front(), map_A);
  auto B = ReadBuchi(input.back(), map_B);
  auto output = DisjointUnion(*A, *B);

  if (verbose != QUIET) {
    output->Print();
  }
}

void print_help(const char *name) {
  printf("Omega Automata Version 6.0 2018/31/30\n");
  printf("Copyright (c) 2011-2018, ");
  printf("Adrian Trejo Nuñez (atrejo@andrew.cmu.edu)\n\n");
  printf("usage: %s [options] file\n", name);
}

int main(int argc, char *argv[]) {
  namespace opt = boost::program_options;

  opt::options_description desc("All Options");
  opt::variables_map var_map;

  desc.add_options()
    ("help,h", "Display This Information")

    ("file,f",
     opt::value<std::vector<std::string>>(),
     "Input File Containing Description of an Automaton")

    ("classify",
     opt::value<std::vector<std::string>>(),
     "Input File Traces of ECAs")

    ("swap,s", "Swap Update and Create in Safra's Algorithm")
    ("clean,c", "Clean Up Rabin Pairs")
    ("minimize,m", "Minimize Rabin Automaton")

    ("product", "Product Construction of Büchi Automata")
    ("union", "Union Construction of Büchi Automata")

    ("canonical",
     opt::value<int>()->notifier(&omega::Canonical),
     "Output Transition Table for Canonical DeBruijn Automaton")
    ("run,r",
     opt::value<int>(),
     "Construct DeBruijn Elementary Cellular Automaton")
    ("cycle",
     opt::value<int>(),
     "Construct a Büchi Automaton to Check for k-Cycles")
    ("dfa",
     opt::value<int>(),
     "Construct Minimal DeBruijn Automaton")
    ("tabulate,t",
     "Output Classification Table of Elementary Cellular Automaton")
    ("k",
     opt::value<int>()->default_value(1),
     "Parameter for Büchi Constructions")
    ("pattern",
     opt::value<std::string>(),
     "Pattern to Avoid in Büchi Constructions")

    ("interactive,i",
     "Interactive Testing")
    ("test,t",
     opt::value<std::string>(),
     "Test Input from File")

    ("outfile,o",
     "Output to .aut File")
    ("dot",
     opt::value<std::string>(),
     "Output To .dot File")
    ("threads",
     opt::value<size_t>(&num_threads)->default_value(4),
     "Number of Threads to Execute Concurrently")
    ("verbose,v",
     opt::value<int>(&verbose)->default_value(GENERAL),
     "Verbosity Level\n  0: Quiet\n  1: General\n  2: Debugging");

  // Allow user to type:
  //   > ./driver file.aut
  // Instead of:
  //   > ./driver --file file.aut
  opt::positional_options_description p;
  p.add("file", -1);

  opt::store(
      opt::command_line_parser(argc, argv).options(desc).positional(p).run(),
      var_map);

  if (var_map.count("help")) {
    print_help(argv[0]);
    std::cout << desc << std::endl;

    std::exit(EXIT_SUCCESS);
  }

  opt::notify(var_map);

  // try {
  //   opt::store(
  //       opt::parse_config_file<char>("omega.cfg", desc), var_map);
  // } catch (const opt::reading_file& e) {
  //   std::cout << "Error: " << e.what() << std::endl;
  // }

  if (var_map.count("dot")) {
    dotfile = var_map["dot"].as<std::string>();
    std::cout << dotfile << std::endl;
  }

  if (var_map.count("outfile")) {
    verbose = OUTFILE;
  }

  std::vector<std::string> patterns;
  if (var_map.count("pattern")) {
    boost::split(
        patterns, var_map["pattern"].as<std::string>(), boost::is_any_of(","));
  }

  // Buchi transition function may be non-deterministic, so allow for multiple
  // transitions under the same symbol.
  TransitionMap map;

  std::unique_ptr<BuchiAutomaton> B;

  std::string filename;
  if (var_map.count("file")) {
    filename = var_map["file"].as<std::vector<std::string>>().front();
    B = ReadBuchi(filename, map);
  } else {
    if (var_map.count("cycle")) {
      if (var_map["k"].as<int>() == 1) {
        omega::FixedPoint(var_map["cycle"].as<int>(), patterns);
      } else {
        omega::Cycle(var_map["cycle"].as<int>(), var_map["k"].as<int>());
      }
    } else if (var_map.count("run")) {
      omega::Run(var_map["run"].as<int>(), var_map["k"].as<int>(), patterns);
    } else if (var_map.count("dfa")) {
      omega::Minimal(var_map["dfa"].as<int>(), var_map["k"].as<int>());
    } else if (var_map.count("tabulate")) {
      omega::Tabulate(var_map["k"].as<int>());
    } else if (var_map.count("product")) {
      Intersection(var_map["file"].as<std::vector<std::string>>());
    } else if (var_map.count("union")) {
      DisjointUnion(var_map["file"].as<std::vector<std::string>>());
    } else if (var_map.count("classify")) {
      filename = var_map["classify"].as<std::vector<std::string>>().front();
      classify(filename);
    } else {
      print_help(argv[0]);
      printf("must specify input file\n");
      std::cout << desc << std::endl;

      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }

  std::unique_ptr<RabinAutomaton> R;
  {
    omega::DeterminizeOpts opts;
    opts.SwapUpdateCreate = var_map.count("swap");
    R = B->Determinize(map, opts);
  }

  if (verbose > QUIET) {
    R->Print();
  }

  if (var_map.count("clean")) {
    R->Clean();
  }

  if (var_map.count("minimize")) {
    R->Minimize();
  }

  assert(!R->pairs.empty());

  if (var_map.count("interactive")) {
    printf("Interactive mode.\n");
    printf("Enter \"u:v\" pair (EOF to stop): ");

    std::string U, V;
    while (ReadUV(std::cin, R->num_alphabet, U, V)) {
      // if input is malformed, ReadUV sets first character of V to nullptr
      if (V[0] != '\0') {
        R->TestUV(U, V);
      }

      printf("\nEnter \"u,v\" pair (EOF to stop): ");
    }

    std::cout << std::endl;
  }

  if (var_map.count("test")) {
    filename = var_map["test"].as<std::string>();

    std::ifstream file(filename);
    if (file.fail()) {
      std::cerr << "Could not open file: " << filename << std::endl;
      exit(EXIT_FAILURE);
    }

    std::string U, V;
    while (ReadUV(file, R->num_alphabet, U, V)) {
      // if input is malformed, ReadUV sets first character of V to nullptr
      if (V[0] != '\0') {
        R->TestUV(U, V);
      }

      std::cout << std::endl;
    }

    file.close();
  }

  return EXIT_SUCCESS;
}
