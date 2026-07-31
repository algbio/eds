CXX_FLAGS=-std=c++17 -O3
#CXX_FLAGS=-std=c++17 -O0 -g
.PHONY : clean
VERSION=$(shell git rev-parse --short HEAD)
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
HTSLIB_INCLUDE=$(ROOT_DIR)ext/htslib
HTSLIB_LIB=$(ROOT_DIR)ext/htslib
HTSLIB_FLAGS=-L $(HTSLIB_LIB) -lhts -Wl,-rpath $(HTSLIB_LIB)
OTHER_INCLUDE=ext/
SDSL_INCLUDE=ext/sdsl-lite/include/

mincard: src/mincard.cpp src/segment.hpp src/rmqueue.h src/rmqueue.cpp src/RMaxQTree.h src/RMaxQTree.cpp src/msa_chunker.hpp src/trie.hpp src/pbwt.h src/pbwt.cpp src/algo.hpp src/minsize.hpp
	${CXX} $(CXX_FLAGS) -DVERSION="\"$(VERSION)\"" -o mincard src/mincard.cpp src/rmqueue.cpp src/RMaxQTree.cpp src/pbwt.cpp -I $(HTSLIB_INCLUDE) -I $(OTHER_INCLUDE) -I $(SDSL_INCLUDE) $(HTSLIB_FLAGS)

clean:
	rm -f mincard

test_rmqueue: src/test_rmqueue.cpp  src/rmqueue.cpp
	${CXX} $(CXX_FLAGS) -o test_rmqueue src/test_rmqueue.cpp src/rmqueue.cpp
