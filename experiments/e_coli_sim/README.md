# `e_coli_sim`
To obtain the dataset, you need [NCBI datasets](https://github.com/ncbi/datasets) and [`iqtree3`](https://iqtree.github.io/). Get them through conda with commands
```sh
conda create -y --name iqtree
conda install -y --name iqtree conda-forge::ncbi-datasets-cli bioconda::iqtree
conda activate iqtree
```
Then, download the `NC_000913.3` E. Coli reference and simulate the MSA with
```sh
./get_dataset.sh
```

Compile `mincard`, obtain the `msatoeds` script, and run the experiment (<100GB of RAM) with commands
```
make -C ../../
git submodule update --init ../ext/junctions
./run_experiment.sh
```
