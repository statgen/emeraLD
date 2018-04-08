VERSION=0.1

CXXFLAGS=-std=c++11 -flto -O4 -pipe

LDFLAGS= -lz -lm

cppsrc = $(wildcard src/*.cpp) src/tabix_util/tabix.cpp
csrc = src/tabix_util/index.c src/tabix_util/bgzf.c

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o)

bin/emeraLD: $(objs)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
