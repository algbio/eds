#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd "$thisfolder"

mincard="$thisfolder/../../mincard"
inputmsa="$thisfolder/input/chr19_100.aligned.uppercase.fa"
usrbintimeformat="%e"

outdir="pbwt_output"
mkdir -p "$outdir"
rm -f "$outdir"/*
cd "$outdir"
ln -s "$inputmsa" msa.fa
ln -s "$thisfolder/input/msa.fa.fai" msa.fa.fai

U_values=(4 8 16)
printf "Running mincard with both suffix trie and pBWT on different U values %s \n" "${U_values[*]}"

tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

declare -a trie_times
declare -a pbwt_times

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m'

# mincard with gaps as symbols strategy
# suffix trie vs pbwt
for U in "${U_values[@]}"
do
    printf "\n${BLUE}Comparison for U = $U${NC}\n"
	# trie run
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard" msa.fa -v --gaps-as-symbols -U $U -o mincard_U${U}.eds
	t1=$(<"$tmpfile")
	trie_times+=("$t1")
	# pbwt run
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard" msa.fa --gaps-as-symbols --pbwt -U $U -o mincard_U${U}_pbwt.eds
	t2=$(<"$tmpfile")
	pbwt_times+=("$t2")

	if cmp -s mincard_U${U}.eds mincard_U${U}_pbwt.eds;
	then
		printf "${GREEN}Outputs are identical for U = ${U}!${NC}\n"
		printf "${PURPLE}Suffix trie time %.3fs vs. pBWT time %.3fs${NC}\n" "$t1" "$t2"
	else 
		printf "${RED}Outputs are different for U = ${U}!${NC}\n"
		exit 1
	fi
done
# Write results to table
echo -e "U value,Trie,pBWT" > results.csv
for i in "${!U_values[@]}";
do
	echo -e "${U_values[$i]},${trie_times[$i]},${pbwt_times[$i]}" >> results.csv
done
printf "\n${BLUE}Wrote algorithm times to results.csv${NC}\n"