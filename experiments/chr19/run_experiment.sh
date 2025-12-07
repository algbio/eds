#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../msa2eds-mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsa=$thisfolder/input/chr19_100.aligned.uppercase.fa
usrbintimeformat="%e total time"
timeouttime="24h"

mkdir output
cd output
ln -s $inputmsa msa.fa

# mincard
for U in 4 8 16 32 64
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa $U
	mv msa.fa.eds mincard_U${U}.eds
done

# mincard pc
for U in 4 8 16 32 64
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa $U 1
	mv msa.fa.eds mincard_U${U}_perfectcols.eds
done

# msatoeds heuristics
for strat in trivial greedy double-greedy
do
	echo "Strategy ${strat}"
	/usr/bin/time -f"$usrbintimeformat" timeout $timeouttime python3 $seqtoed msa.fa "${strat}.eds" ${strat} || true
	if [ -e "eds_${strat}.txt" ]
	then
		python3 $getstats "${strat}.eds" eds
	else
		echo "Cannot compute EDS stats, file not found (possibly due to timeout)"
	fi
done
