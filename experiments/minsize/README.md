# min range experiment

Simple minsize experiment to compare the different ways of calculating the range minimum query for the M matrix rows.

Get covid dataset as described [here](../covid/README.md)

Compile mincard with each algorithm (ring buffer, rmaxtree, rmqueue) by changing [this line](../../src/minsize.hpp#L25) and then running

```sh
 make -C ../../
 mv ../../mincard ./mincard-ring
```

Do this with the names [ring, tree, queue]. Finally run the experiment script (running times will be collected in output/results.csv)

```sh
 ./run_experiment.sh
```
