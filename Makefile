VERSION=0.1

CPPFLAGS= -Ofast -flto -pipe -Isrc
CXXFLAGS= -std=c++11

LDFLAGS= -lz -lm

.PHONY: all clean

cppsrc = $(wildcard src/*.cpp) src/tabix_util/tabix.cpp
csrc = src/tabix_util/index.c src/tabix_util/bgzf.c

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

