# max range experiment

Simple pbwt experiment to compare the different ways of calculating the range maximum query for the ek array.

Get covid dataset as described [here](../covid/README.md)

Compile mincard with each algorithm (naive, recursive, rmq) by changing [this line](../../src/pbwt.hpp#L35) and then running

```sh
 make -C ../../
 mv ../../mincard ./mincard-naive
```

Finally run the experiment script (running times will be collected in output/results.csv)

```sh
 ./run_experiment.sh
```
