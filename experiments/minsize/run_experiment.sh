#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd "$thisfolder"

mincard_ring="$thisfolder/mincard-ring"
mincard_tree="$thisfolder/mincard-tree"
mincard_queue="$thisfolder/mincard-queue"
inputmsa="$thisfolder/input/covid_100000.fa.gz"
usrbintimeformat="%e"

outdir="output"
mkdir -p "$outdir"
rm -f "$outdir"/*
cd "$outdir"
ln -s "$inputmsa" msa.fa

L_values=(1 2 4 8 16 32 64 128 256 512)
printf "Running minsize with pBWT and different range min strategies (ring, rmaxtree, rmqueue) on different L values %s \n" "${L_values[*]}"

tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m'

# Write results to table
echo -e "L value,Ring Buffer,RMaxTree,RMQueue" > results.csv

# minsize with pbwt (gaps as symbols strategy)
for L in "${L_values[@]}"
do
	printf "\n${BLUE}Comparison for L = $L${NC}\n"
	# ring buffer
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard_ring" msa.fa -v --gaps-as-symbols --pbwt -L $L --min-size -o minsize_ring_L${L}.eds
	t1=$(<"$tmpfile")
	# rmaxtree
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard_tree" msa.fa -v --gaps-as-symbols --pbwt -L $L --min-size -o minsize_tree_L${L}.eds
	t2=$(<"$tmpfile")
	# rmqueue
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard_queue" msa.fa -v --gaps-as-symbols --pbwt -L $L --min-size -o minsize_queue_L${L}.eds
	t3=$(<"$tmpfile")
	
	if cmp -s minsize_ring_L${L}.eds minsize_tree_L${L}.eds &&
			cmp -s minsize_ring_L${L}.eds minsize_queue_L${L}.eds
	then
		printf "${GREEN}Outputs are identical for L = ${L}!${NC}\n"
		printf "${PURPLE}ring buffer %.3fs vs. rmaxtree suffix trie %.3fs vs. rmqueue %.3fs${NC}\n" "$t1" "$t2" "$t3"
	else 
		printf "${RED}Outputs are different for ${name} on U = ${U}!${NC}\n"
		exit 1
	fi

	echo -e "${L},${t1},${t2},${t3}" >> results.csv
done

printf "\n${GREEN}Finished running the experiment!\n${BLUE}Wrote algorithm times to results.csv${NC}\n"
