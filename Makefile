CXX_FLAGS=-std=c++17 -O3
#CXX_FLAGS=-std=c++17 -O0 -g
.PHONY : clean
VERSION=$(shell git rev-parse --short HEAD)
HTSLIB_INCLUDE=ext/htslib
HTSLIB_LIB=ext/htslib
HTSLIB_FLAGS=-L $(HTSLIB_LIB) -lhts -Wl,-rpath,$(HTSLIB_LIB)

msa2eds-mincard: src/msa2eds-mincard.cpp src/block_graph.hpp src/RMaxQTree.h src/RMaxQTree.cpp src/msa_chunker.hpp
	${CXX} $(CXX_FLAGS) -DVERSION="\"$(VERSION)\"" -o msa2eds-mincard src/msa2eds-mincard.cpp src/RMaxQTree.cpp -I $(HTSLIB_INCLUDE) $(HTSLIB_FLAGS)

clean:
	rm -f msa2eds-mincard
