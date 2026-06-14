# Build configuration: debug (default) or release
#
#   make BUILD=release
#
# Debug builds are unoptimized, carry debug info, and define DEBUG_BUILD to
# enable the dbg() logging macros
# Release builds are optimized and define NDEBUG
BUILD ?= debug

ifeq ($(BUILD),release)
  OPTFLAGS = -O2
  BUILD_DEFS = -DNDEBUG
else ifeq ($(BUILD),debug)
  OPTFLAGS = -O0 -g
  BUILD_DEFS = -DDEBUG_BUILD
else
  $(error Unknown BUILD '$(BUILD)'; use 'debug' or 'release')
endif

# Boost location
#
# Override any of these on the command line, e.g.
#   make BOOST_PREFIX=/usr/local
#   make BOOST_INC=/path/to/include BOOST_LIB=/path/to/lib
#
# By default we try to auto-detect the prefix: the Homebrew prefix on macOS
# (via brew, falling back to /opt/homebrew) and the standard system prefix on
# Linux, where libboost-all-dev installs into /usr
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	BOOST_PREFIX ?= $(shell brew --prefix 2>/dev/null || (test -d /opt/homebrew && echo /opt/homebrew) || echo /usr/local)
else
	BOOST_PREFIX ?= /usr
endif

BOOST_INC ?= $(BOOST_PREFIX)/include
BOOST_LIB ?= $(BOOST_PREFIX)/lib

# Only emit -I/-L for non-standard prefixes so system installs work cleanly
BOOST_CXXFLAGS = $(if $(filter-out /usr/include,$(BOOST_INC)),-isystem $(BOOST_INC))
BOOST_LDFLAGS = $(if $(filter-out /usr/lib,$(BOOST_LIB)),-L$(BOOST_LIB))

CC = clang
CFLAGS = -Wall -Wextra $(OPTFLAGS) $(BUILD_DEFS)

CXX = clang++
CXXFLAGS = --std=c++23 -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter \
					$(OPTFLAGS) $(BUILD_DEFS) \
					-I. $(BOOST_CXXFLAGS)

					# -Wthread-safety \
					# -Wfor-loop-analysis \
					# -fno-elide-type \
					# -fdiagnostics-show-template-tree
					# -Weverything \
					# -Wno-c++98-compat -Wno-c++98-compat-pedantic \
					# -Wno-c++11-compat -Wno-c++11-compat-pedantic \
					# -Wno-c++14-compat -Wno-c++14-compat-pedantic \
					# -Wno-missing-prototypes \
					# -Wno-missing-noreturn \
					# -Wno-format-nonliteral \

# Linker search paths / flags only
# Library flags go in LDLIBS so they can be placed after the object files
# on the link line: GNU ld with --as-needed
# (the default on many Linux toolchains) drops libraries that appear before
# the objects referencing them, causing undefined-reference errors
LDFLAGS = $(BOOST_LDFLAGS)

LDLIBS = -lboost_program_options -lboost_stacktrace_basic

					# -fsanitize=thread,undefined,integer,nullability,safe-stack \
					# -fsanitize-recover=all
					# -ltcmalloc -lprofiler -Wl,-no_pie

.SUFFIXES:
.SUFFIXES: .cpp .h .o

SHELL = /bin/sh
RM = rm -f

# Per-configuration directory for object files and other intermediate build
# artifacts, so that switching BUILD never reuses objects compiled with a
# different configuration
# The final binaries are written to the repo root
BUILDDIR = build/$(BUILD)

# Auto-generate header dependencies (.d files) alongside each object so that
# changing a header only rebuilds the objects that actually include it
DEPFLAGS = -MMD -MP

SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst %.cpp,$(BUILDDIR)/%.o,$(notdir $(SRCS)))
TGTS = driver test

# Objects shared by every executable: every compiled source except the
# per-executable entry points
# Derive these from the source list so a header-only file never makes Make
# chase a non-existent object
EXEC_OBJS = $(patsubst %,$(BUILDDIR)/%.o,$(TGTS))
LIB_OBJS = $(filter-out $(EXEC_OBJS),$(OBJS))

.PHONY: all clean log

all: $(TGTS)

$(BUILDDIR):
	mkdir -p $@

batch: src/batch.c
	$(CC) $(CFLAGS) $< -o $@

$(BUILDDIR)/%.o : src/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

driver: $(LIB_OBJS) $(BUILDDIR)/driver.o
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

test: $(LIB_OBJS) $(BUILDDIR)/test.o
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

log:
	VERBOSE=3 ./test >& log

clean:
	$(RM) -r build $(TGTS) batch

-include $(OBJS:.o=.d)
