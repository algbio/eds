#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../msa2eds-mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsa=$thisfolder/input/covid19-100-N.fa
usrbintimeformat="%e total time"

mkdir output
cd output
ln -s $inputmsa msa.fa

# mincard
for U in 4 8 16 32 64 128 256 512
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa $U
	mv msa.fa.eds mincard_U${U}.eds
done

# mincard pc
for U in 4 8 16 32 64 128 256 512
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa $U 1
	mv msa.fa.eds mincard_U${U}_perfectcols.eds
done

# mincard trivial S^|||
/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa 0 0 1

# msatoeds heuristics
for strat in trivial greedy double-greedy
do
	echo "Strategy ${strat}"
	/usr/bin/time -f"$usrbintimeformat" python3 $seqtoed msa.fa "${strat}.eds" ${strat}
	python3 $getstats "${strat}.eds" eds
done
