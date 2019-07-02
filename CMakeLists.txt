cmake_minimum_required(VERSION 3.9)
project(emeraLD LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)

find_library(LIB_Z z)
find_library(LIB_M m)

include_directories(src src/tabix_util)

add_library(tabix_util
        src/tabix_util/bam_endian.h
        src/tabix_util/bedidx.c
        src/tabix_util/bgzf.c
        src/tabix_util/bgzf.h
        src/tabix_util/index.c
        src/tabix_util/khash.h
        src/tabix_util/knetfile.c
        src/tabix_util/knetfile.h
        src/tabix_util/kseq.h
        src/tabix_util/ksort.h
        src/tabix_util/kstring.c
        src/tabix_util/kstring.h
        src/tabix_util/tabix.cpp
        src/tabix_util/tabix.h
        src/tabix_util/tabix.hpp
)

set_target_properties(tabix_util PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
target_include_directories(tabix_util PUBLIC src/tabix_util)
target_link_libraries(tabix_util ${LIB_Z} ${LIB_M})

add_executable(emeraLD
        src/calcLD.cpp
        src/calcLD.hpp
        src/Main.cpp
        src/processGenotypes.cpp
        src/processGenotypes.hpp)

target_link_libraries(emeraLD tabix_util)

add_executable(emeraldpy
        src/emeraldpy/core.cpp
        src/calcLD.cpp
        src/processGenotypes.cpp
)

target_link_libraries(emeraldpy tabix_util)

# Specifically needed on MacOS.
#set_target_properties(emeraldpy PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")