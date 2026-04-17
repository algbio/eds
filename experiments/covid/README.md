# `covid` experiment
To obtain the datasets, you need [NCBI datasets](https://github.com/ncbi/datasets), [`minimap2`](https://github.com/lh3/minimap2), and [`ViralMSA`](https://github.com/niemasd/ViralMSA). Get `dataset` through conda and the other dependencies with commands
```sh
conda create -y --name ncbi
conda install -y --name ncbi -c conda-forge ncbi-datasets-cli
conda activate ncbi
git submodule update --init --recursive ../ext/minimap2 ../ext/ViralMSA
make -C ../ext/minimap2
```
Finally, the sequences can be downloaded and aligned with
```sh
./get_datasets.sh
```

Compile `mincard`, obtain the `msatoeds` script, and run the experiment with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```

Afterwards, the resulting EDSes can be (slowly) verified using 16 threads with script
```
./verify_edses.sh 16
```
