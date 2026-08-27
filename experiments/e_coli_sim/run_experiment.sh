#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsa=$thisfolder/input/msa.fa
usrbintimeformat="%e total time"

mkdir output
cd output
ln -s $inputmsa msa.fa

# mincard
for U in 4 8 16 32 64
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa -v --gaps -U $U -o mincard_U${U}.eds
done

# mincard with perfect segments
for U in 4 8 16 32 64
do
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa -v --gaps -U $U --perfect-segments -o mincard_U${U}_p.eds
done

# mincard trivial S^|||
/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa -v --gaps --trivial-vertical -o mincard_t.eds

# mincard trivial S^≡ with perfect segments
/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa -v --gaps --trivial-horizontal --perfect-segments -o mincard_np.eds

# msatoeds heuristics
for strat in trivial greedy double-greedy
do
	echo "Strategy ${strat}"
	/usr/bin/time -f"$usrbintimeformat" python3 $seqtoed msa.fa "${strat}.eds" ${strat}
	python3 $getstats "${strat}.eds" eds
done
