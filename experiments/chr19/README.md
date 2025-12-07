# `chr19`
After checking out the "getting the dataset" subsection, compile `msa2eds-mincard`, obtain the `msatoeds` scripts, and run the experiment with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```

Afterwards, the resulting EDSes can be verified using 16 threads with script
```
./verify_edses.sh 16
```

## getting the dataset
The script `get_datasets.sh` assumes that you have [`seqtk`](https://github.com/lh3/seqtk) installed and visible by your `PATH` environment variable (modify `get_datasets.sh` accordingly if you get the program some other way).
```
./get_datasets.sh
```

