#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../msa2eds-mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsa=$thisfolder/input/msa.fa
usrbintimeformat="%e total time"

mkdir output
cd output
ln -s $inputmsa msa.fa

# mincard
for U in 4 8 16 32 64 128
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa $U
	mv msa.fa.gfa eds_U$U.gfa
done

# mincard pc
for U in 4 8 16 32 64 128
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa $U 1
	mv msa.fa.gfa eds_U$U.gfa
done

# msatoeds heuristics
for strat in trivial greedy double-greedy
do
	echo "Strategy ${strat}"
	/usr/bin/time -f"$usrbintimeformat" python3 $seqtoed msa.fa "eds_${strat}.txt" ${strat}
	python3 $getstats "eds_${strat}.txt" eds
done
