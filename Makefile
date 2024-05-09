VERSION = 5.0

CC = clang
CFLAGS = -Wall -Wextra

CPPFLAGS = -DDEBUG_BUILD
CXX = clang++
CXXFLAGS = --std=c++17 -Wall -Wextra \
					 -I. -isystem /usr/local/include \
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

LDFLAGS = -L/usr/local/lib \
					-lboost_program_options-mt \
					-lboost_stacktrace_basic-mt

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
