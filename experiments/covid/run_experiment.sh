#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsas=($thisfolder/input/covid_100.fa.gz $thisfolder/input/covid_1000.fa.gz $thisfolder/input/covid_10000.fa.gz $thisfolder/input/covid_100000.fa.gz)
format="LOG: %e total time (s), %M maximum resident set size (KB)\n"

mkdir output
cd output

for inputmsa in "${inputmsas[@]}"
do
	if [ ! -f $inputmsa ]
	then
		echo "ERROR: input dataset $inputmsa not found!"; exit 1
	fi
done

for inputmsa in "${inputmsas[@]}"
do
	base=$(basename $inputmsa .fa.gz)
	echo "LOG: processing $base..."
	ln -s $inputmsa msa.fa.gz # create symlink

	for U in 4 8 16 32 64 128 256 512
	do
		# mincard OPT
		/usr/bin/time -f"$format" $mincard \
			msa.fa.gz -o ${base}_mincard_U${U}.eds \
			--no-pbwt --gaps-as-gaps --disable-perfect-segments --max-segment-length $U --verbose
		rm msa.fa.gz.fai msa.fa.gz.gzi

		# mincard OPT pc
		/usr/bin/time -f"$format" $mincard \
			msa.fa.gz -o ${base}_mincard_U${U}_p.eds \
			--no-pbwt --gaps-as-gaps --max-segment-length $U --verbose
		rm msa.fa.gz.fai msa.fa.gz.gzi

		# mincard OPT_-∈Σ (pBWT)
		/usr/bin/time -f"$format" $mincard \
			msa.fa.gz -o ${base}_mincard_U${U}_g.eds \
			--disable-perfect-segments --max-segment-length $U --verbose
		rm msa.fa.gz.fai msa.fa.gz.gzi

		# mincard OPT_-∈Σ pc (pBWT)
		/usr/bin/time -f"$format" $mincard \
			msa.fa.gz -o ${base}_mincard_U${U}_pg.eds \
			--max-segment-length $U --verbose
		rm msa.fa.gz.fai msa.fa.gz.gzi
	done

	for U in 1024 2048 4096 8192
	do
		# mincard OPT_-∈Σ (pBWT)
		/usr/bin/time -f"$format" $mincard \
			msa.fa.gz -o ${base}_mincard_U${U}_g.eds \
			--disable-perfect-segments --max-segment-length $U --verbose
		rm msa.fa.gz.fai msa.fa.gz.gzi

		# mincard OPT_-∈Σ pc (pBWT)
		/usr/bin/time -f"$format" $mincard \
			msa.fa.gz -o ${base}_mincard_U${U}_pg.eds \
			--max-segment-length $U --verbose
		rm msa.fa.gz.fai msa.fa.gz.gzi
	done

	# mincard trivial S^|||
	/usr/bin/time -f"$format" $mincard \
		msa.fa.gz -o ${base}_mincard_t.eds \
		--verbose --trivial-vertical --disable-perfect-segments
	rm msa.fa.gz.fai msa.fa.gz.gzi

	# mincard trivial S^≡ pc
	/usr/bin/time -f"$format" $mincard \
		msa.fa.gz -o ${base}_mincard_np.eds \
		--verbose --trivial-horizontal
	rm msa.fa.gz.fai msa.fa.gz.gzi

	# msatoeds heuristics
	for strat in trivial greedy
	do
		echo "Strategy ${strat}"
		/usr/bin/time -f"$format" python3 $seqtoed <(gunzip -c msa.fa.gz) "${base}_${strat}.eds" ${strat}
		python3 $getstats "${base}_${strat}.eds" eds
	done

	rm msa.fa.gz # remove symlink
done
