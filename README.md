# eds
Program `msa2eds-mincard` constructs Elastic Degenerate Strings from a multiple sequence alignment in FASTA format. Tested on GCC 15.2.1. You can get this repository and compile the [HTSlib](https://github.com/samtools/htslib) dependency with
```sh
git clone https://github.com/algbio/eds && cd eds
git submodule update --init --recursive ext/htslib
(cd ext/htslib && autoreconf -i && ./configure && make -j)
```
Otherwise, if you install `htslib` some other way (for example from its [website](https://www.htslib.org/download/)), change variables `X` and `Y` in `Makefile`. Finally, compile `msa2esd-mincard` and test it with
```sh
make
./msa2eds-mincard test/example.fasta
```

## todo
- add option to preprocess meaningful extensions
- add verbose option logging times
- fix tests
- support bgzipped files
- QC on the edses in chr19 experiment
