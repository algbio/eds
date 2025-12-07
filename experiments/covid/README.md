# `covid`
Get the dataset with command
```
./get_datasets.sh
```

Compile `msa2eds-mincard`, obtain the `msatoeds` scripts, and run the experiment with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```

Afterwards, the resulting EDSes can be verified using 16 threads with script
```
./verify_edses.sh 16
```
