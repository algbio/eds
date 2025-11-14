FLAGS=-std=c++17 -O2
.PHONY : clean

msa2eds-mincard:
	${CXX} $(FLAGS) src/msa2eds-mincard.cpp src/RMaxQTree.cpp -o msa2eds-mincard

clean:
	rm -f msa2eds-mincard
