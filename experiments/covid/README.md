# `covid` experiment
To obtain the datasets, you need [NCBI datasets](https://github.com/ncbi/datasets), [`minimap2`](https://github.com/lh3/minimap2), and [`ViralMSA`](https://github.com/niemasd/ViralMSA). Get `dataset` through conda and the other dependencies with commands
```sh
conda create -y --name ncbi
conda install -y --name ncbi -c conda-forge ncbi-datasets-cli bioconda::bcftools
conda activate ncbi
git submodule update --init --recursive ../ext/minimap2 ../ext/ViralMSA
make -C ../ext/minimap2
```
Finally, the sequences can be downloaded and aligned (~4GB of disk space and <16GB of RAM) with
```sh
./get_dataset.sh
```

Compile `mincard`, obtain the `msatoeds` script, and run the experiment (<100GB of RAM) with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```
