# `chr19`
After checking out the "getting the dataset" subsection, compile `msa2eds-mincard`, obtain the `msatoeds` scripts, and run the experiment with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```

## getting the dataset
The script `get_datasets.sh` assumes that you have [`seqtk`](https://github.com/lh3/seqtk) installed and visible by your `PATH` environment variable (modify `get_datasets.sh` accordingly if you get the programs some other way).
```
./get_datasets.sh
```

