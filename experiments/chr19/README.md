# `chr19` experiment
To obtain the dataset, you need [`seqtk`](https://github.com/lh3/seqtk). If it's not installed on your system, get it and compile it with commands
```sh
git submodule update --init ../ext/seqtk
make -C ../ext/seqtk
```
Then, get the dataset (~60GB) with command
```sh
./get_dataset.sh
```

Compile `mincard`, obtain the `msatoeds` script, and run the experiment (<100GB of RAM) with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```
