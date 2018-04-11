VERSION=0.1

CPPFLAGS=-Ofast -pipe -Wall
CXXFLAGS=-std=c++11

LDFLAGS= -lz -lm

cppsrc = $(wildcard src/*.cpp) src/tabix_util/tabix.cpp
csrc = src/tabix_util/index.c src/tabix_util/bgzf.c

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o)

bin/emeraLD: $(objs)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean: 
	rm -f src/*.o src/*/*.o bin/emeraLD

