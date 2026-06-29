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

mincard: src/mincard.cpp src/segment.hpp src/RMaxQTree.h src/RMaxQTree.cpp src/msa_chunker.hpp src/trie.hpp src/pbwt.hpp src/algo.hpp
	${CXX} $(CXX_FLAGS) -DVERSION="\"$(VERSION)\"" -o mincard src/mincard.cpp src/RMaxQTree.cpp -I $(HTSLIB_INCLUDE) -I $(OTHER_INCLUDE) -I $(SDSL_INCLUDE) $(HTSLIB_FLAGS)

clean:
	rm -f mincard
