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

msa2eds-mincard: src/msa2eds-mincard.cpp src/block_graph.hpp src/RMaxQTree.h src/RMaxQTree.cpp src/msa_chunker.hpp src/trie.hpp
	${CXX} $(CXX_FLAGS) -DVERSION="\"$(VERSION)\"" -o msa2eds-mincard src/msa2eds-mincard.cpp src/RMaxQTree.cpp -I $(HTSLIB_INCLUDE) -I $(OTHER_INCLUDE) -I $(SDSL_INCLUDE) $(HTSLIB_FLAGS)

clean:
	rm -f msa2eds-mincard
