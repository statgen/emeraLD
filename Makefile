VERSION=0.1

CPPFLAGS= -Ofast -flto -pipe -Isrc
CXXFLAGS= -std=c++11

LDFLAGS= -lz -lm

.PHONY: all clean

cppsrc = $(wildcard src/*.cpp) src/tabixpp/tabix.cpp
csrc = src/htslib/index.c src/htslib/bgzf.c

all: version

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o)

version:
	@echo "Building version $(VERSION)"
	@echo "#define VERSION $(VERSION)" > src/version.hpp
	@$(MAKE) --no-print-directory bin/emeraLD

bin/emeraLD: $(objs)
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

clean: 
	rm -f src/*.o src/*/*.o bin/emeraLD

