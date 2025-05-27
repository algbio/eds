# eds
Constructing EDSes from MSAs
# compile 
g++ -std=c++17 -O2 msa2eds-mincard.cpp RMaxQTree.cpp -o msa2eds-mincard
# execute
./msa2eds-mincard example.fasta 4
