VERSION = 5.0

CC = clang
CFLAGS = -Wall -Wextra

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
BOOST_CPPFLAGS = $(if $(filter-out /usr/include,$(BOOST_INC)),-isystem $(BOOST_INC))
BOOST_LDFLAGS = $(if $(filter-out /usr/lib,$(BOOST_LIB)),-L$(BOOST_LIB))

CPPFLAGS = -DDEBUG_BUILD
CXX = clang++
CXXFLAGS = --std=c++23 -Wall -Wextra \
					 -I. $(BOOST_CPPFLAGS) \
					 -Wno-unused-variable -Wno-unused-parameter -g

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

LDFLAGS = $(BOOST_LDFLAGS) \
					-lboost_program_options \
					-lboost_stacktrace_basic

					# -fsanitize=thread,undefined,integer,nullability,safe-stack \
					# -fsanitize-recover=all
					# -ltcmalloc -lprofiler -Wl,-no_pie

.SUFFIXES:
.SUFFIXES: .cpp .h .o

SHELL = /bin/sh
RM = rm -f

TMPDIR = build

SRCS = $(wildcard src/*.cpp)
HDRS = $(wildcard src/*.h)
DEPS = $(patsubst %.h,$(TMPDIR)/%.o,$(notdir $(HDRS)))
OBJS = $(patsubst %.cpp,$(TMPDIR)/%.o,$(notdir $(SRCS)))
TGTS = driver test

$(info $$SRCS = ${SRCS})
$(info $$HRDS = ${HDRS})
$(info $$DEPS = ${DEPS})
$(info $$OBJS = ${OBJS})
$(info $$TGTS = ${TGTS})
$(shell mkdir -p $(TMPDIR))

.PHONY: all clean log

all: $(TGTS)

batch:
	$(CC) $(CFLAGS) batch.c -o $(TMPDIR)/$@

$(TMPDIR)/%.o : src/%.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

driver: $(DEPS) $(TMPDIR)/driver.o
	$(CXX) $(LDFLAGS) $^ -o $@

test: $(DEPS) $(TMPDIR)/test.o
	$(CXX) $(LDFLAGS) $^ -o $@

log:
	VERBOSE=3 ./test >& log

clean:
	$(RM) $(OBJS) $(TGTS)
