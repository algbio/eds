# eds
Program `mincard` constructs Elastic Degenerate Strings from a multiple sequence alignment in FASTA format. Tested on GCC 15.2.1. You can get this repository and compile the [HTSlib](https://github.com/samtools/htslib) and [SDSL v3](https://github.com/xxsds/sdsl-lite) dependencies with
```sh
git clone https://github.com/algbio/eds && cd eds
git submodule update --init --recursive ext/htslib ext/sdsl-lite
(cd ext/htslib && autoreconf -i && ./configure && make)
```
Otherwise, if you install HTSlib some other way (for example from its [website](https://www.htslib.org/download/)), change variables `HTSLIB_INCLUDE` and `HTSLIB_LIB` in `Makefile`. Finally, compile `mincard` and test it with
```sh
make
./mincard test/msa.fa -o test/msa.eds
```

## todo
- compact tries to reduce topology ops
- QC on the EDSes in experiment (maybe with https://github.com/giovannarosone/EDS-BWT?)
- show VCF workflow in README
