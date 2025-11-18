# `e_coli_sim`
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
