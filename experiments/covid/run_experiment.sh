#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsas=($thisfolder/input/covid_100.fa.gz $thisfolder/input/covid_1000.fa.gz $thisfolder/input/covid_10000.fa.gz $thisfolder/input/covid_100000.fa.gz)
usrbintimeformat="%e total time"

mkdir output
cd output

for inputmsa in "${inputmsas[@]}"
do
	base=$(basename $inputmsa .fa.gz)
	echo "LOG: processing $base..."
	ln -s $inputmsa msa.fa.gz
	# mincard
	for U in 4 8 16 32 64 128 256 512
	do
		/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa.gz -v --gaps -U $U -o ${base}_mincard_U${U}.eds # plain
		rm msa.fa.gz.fai msa.fa.gz.gzi
		/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa.gz -v --gaps -U $U --perfect-segments -o ${base}_mincard_U${U}_p.eds # perfect segments
		rm msa.fa.gz.fai msa.fa.gz.gzi
		/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa.gz --gaps --preprocess -v -U $U -o ${base}_mincard_U${U}.eds # preprocess
		rm msa.fa.gz.fai msa.fa.gz.gzi
		/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa.gz -v --gaps -U $U --perfect-segments --preprocess -o ${base}_mincard_U${U}_p.eds
		rm msa.fa.gz.fai msa.fa.gz.gzi
	done

	# mincard trivial S^|||
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa.gz -v --gaps --trivial-vertical -o ${base}_mincard_t.eds

	# mincard trivial S^≡ with perfect segments
	/usr/bin/time -f"$usrbintimeformat" $mincard msa.fa.gz -v --gaps --trivial-horizontal --perfect-segments -o ${base}_mincard_np.eds

	# msatoeds heuristics
	for strat in trivial greedy double-greedy
	do
		echo "Strategy ${strat}"
		/usr/bin/time -f"$usrbintimeformat" python3 $seqtoed <(gunzip -c msa.fa.gz) "${base}_${strat}.eds" ${strat}
		python3 $getstats "${base}_${strat}.eds" eds
	done
	rm msa.fa.gz
done
