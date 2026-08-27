# cardinality vs size experiment

This experiment compares the resulting eds stats and running times for all 4 algorithms (mincard/minsize x gaps/gapless) with various bounds (L and U) on all 3 datasets.

Compile `mincard` and `minsize`:
```
make -C ../../
make minsize -C ../../
```

Get the ecoli, covid, and chr datasets in their respective experiments. Then, simply run the experiment. 

This experiment does not stream the eds-es, though it could with a minor change in the script. The results are collected automatically in output/results.csv
