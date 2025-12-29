FLAGS=-std=c++17 -O3
#FLAGS=-std=c++17 -O0 -g
.PHONY : clean
VERSION=$(shell git rev-parse --short HEAD)

msa2eds-mincard: src/msa2eds-mincard.cpp src/block_graph.hpp src/RMaxQTree.h src/RMaxQTree.cpp
	${CXX} $(FLAGS) src/msa2eds-mincard.cpp src/RMaxQTree.cpp -DVERSION="\"$(VERSION)\"" -o msa2eds-mincard

clean:
	rm -f msa2eds-mincard
