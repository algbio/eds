# column-major experiment

This experiment compares the running times of the suffix tree and pbwt algorithms on both row-major and column-major large msas.

Get covid and chr datasets in their respective experiments. Then, transpose the msas and run `./run_experiment.sh`

Transpose preprocessing times:

- covid_10000: 2.55s
- covid_100000: 43.71s
- chr19: 765.18s

The results indicate that processing the msa in a column-major fashion only benefits the pBWT algorithm and there's an improvement only when the number of rows is large.
