#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../mincard
seqtoed=$thisfolder/../ext/junctions/scripts/msatoeds/seq_to_ed.py
getstats=$thisfolder/../ext/junctions/scripts/msatoeds/get_stats.py
inputmsa=$thisfolder/input/msa.fa
format="LOG: %e total time (s), %M maximum resident set size (KB)\n"

mkdir output
cd output
ln -s $inputmsa msa.fa # create symlink

for U in 4 8 16 32 64
do
	# mincard OPT
	/usr/bin/time -f"$format" $mincard \
		msa.fa -o mincard_U${U}.eds \
		--no-pbwt --gaps-as-gaps --disable-perfect-segments --max-segment-length $U --verbose
	rm msa.fa.fai

	# mincard OPT pc
	/usr/bin/time -f"$format" $mincard \
		msa.fa -o mincard_U${U}_p.eds \
		--no-pbwt --gaps-as-gaps --max-segment-length $U --verbose
	rm msa.fa.fai

	# mincard OPT_-∈Σ (pBWT)
	/usr/bin/time -f"$format" $mincard \
		msa.fa -o mincard_U${U}_g.eds \
		--disable-perfect-segments --max-segment-length $U --verbose
	rm msa.fa.fai

	# mincard OPT_-∈Σ pc (pBWT)
	/usr/bin/time -f"$format" $mincard \
		msa.fa -o mincard_U${U}_pg.eds \
		--max-segment-length $U --verbose
	rm msa.fa.fai
done

for U in 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288
do
	# mincard OPT_-∈Σ (pBWT)
	/usr/bin/time -f"$format" $mincard \
		msa.fa -o mincard_U${U}_g.eds \
		--disable-perfect-segments --max-segment-length $U --verbose
	rm msa.fa.fai

	# mincard OPT_-∈Σ pc (pBWT)
	/usr/bin/time -f"$format" $mincard \
		msa.fa -o mincard_U${U}_pg.eds \
		--max-segment-length $U --verbose
	rm msa.fa.fai
done

# mincard trivial S^|||
/usr/bin/time -f"$format" $mincard \
	msa.fa -o mincard_t.eds \
	--verbose --trivial-vertical --disable-perfect-segments
rm msa.fa.fai

# mincard trivial S^≡ pc
/usr/bin/time -f"$format" $mincard \
	msa.fa -o mincard_np.eds \
	--verbose --trivial-horizontal
rm msa.fa.fai

# msatoeds heuristics
for strat in trivial greedy
do
	echo "Strategy ${strat}"
	/usr/bin/time -f"$format" python3 $seqtoed msa.fa "${strat}.eds" ${strat}
	python3 $getstats "${strat}.eds" eds
done
