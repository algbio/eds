#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$thisfolder"
export LC_NUMERIC="en_US.UTF-8"

mincard="$thisfolder/../../mincard"

# MSA fasta files and column major txt files
input_covid10000_msa="$thisfolder/../covid/input/covid_10000.fa.gz"
input_covid10000_cm_msa="$thisfolder/../covid/input/covid_10000.txt"
input_covid100000_msa="$thisfolder/../covid/input/covid_100000.fa.gz"
input_covid100000_cm_msa="$thisfolder/../covid/input/covid_100000.txt"
input_chr19_msa="$thisfolder/../chr19/input/chr19_100.aligned.uppercase.fa"
input_chr19_cm_msa="$thisfolder/../chr19/input/chr19_100.aligned.uppercase.txt"

# Create links to msa files
outdir="output"
mkdir -p "$outdir"
rm -f "$outdir"/*
cd "$outdir"

ln -s "$input_covid10000_msa" covid10000_msa.fa
ln -s "$input_covid10000_cm_msa" covid10000_cm_msa.txt
ln -s "$input_covid100000_msa" covid100000_msa.fa
ln -s "$input_covid100000_cm_msa" covid100000_cm_msa.txt
ln -s "$input_chr19_msa" chr19_msa.fa
ln -s "$input_chr19_cm_msa" chr19_cm_msa.txt

# Keep them in an array for easy iteration
declare -A msa
declare -A cm_msa

datasets=(
    covid10000
    covid100000
    chr19
)

msa[covid10000]=covid10000_msa.fa
cm_msa[covid10000]=covid10000_cm_msa.txt
msa[covid100000]=covid100000_msa.fa
cm_msa[covid100000]=covid100000_cm_msa.txt
msa[chr19]=chr19_msa.fa
cm_msa[chr19]=chr19_cm_msa.txt

usrbintimeformat="%e"

U_values=(4 8 16 32 64)
printf "Running mincard with both suffix trie and pBWT on row/column-major MSAs with U values %s \n" "${U_values[*]}"

tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m'

# Write results to table
echo -e "Dataset,U value,Suffix trie,CM suffix trie,pBWT,CM pBWT" > results.csv

# mincard with gaps as symbols strategy
# suffix trie vs pbwt
# row-major vs column-major
for name in "${datasets[@]}"; do
    input_msa_row="${msa[$name]}"
    input_msa_col="${cm_msa[$name]}"

    for U in "${U_values[@]}"; do
	    printf "\n${BLUE}Dataset: ${name} | U = ${U}${NC}\n"

		# suffix trie
		/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard" "$input_msa_row" -v --trie -U $U -o mincard_trie_U${U}.eds
		t1=$(<"$tmpfile")
		/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard" "$input_msa_col" -v --trie --column-major -U $U -o mincard_cm_trie_U${U}.eds
		t2=$(<"$tmpfile")
		# pbwt
		/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard" "$input_msa_row" -v -U $U -o mincard_pbwt_U${U}.eds
		t3=$(<"$tmpfile")
		/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard" "$input_msa_col" -v --column-major -U $U -o mincard_cm_pbwt_U${U}.eds
		t4=$(<"$tmpfile")
	
		if cmp -s mincard_trie_U${U}.eds mincard_cm_trie_U${U}.eds &&
			cmp -s mincard_trie_U${U}.eds mincard_pbwt_U${U}.eds &&
			cmp -s mincard_trie_U${U}.eds mincard_cm_pbwt_U${U}.eds;	
		then
			printf "${GREEN}Outputs are identical for ${name} on U = ${U}!${NC}\n"
			printf "${PURPLE}suffix trie %.3fs vs. column-major suffix trie %.3fs vs. pbwt %.3fs vs. column-major pbwt %.3fs${NC}\n" "$t1" "$t2" "$t3" "$t4"
		else 
			printf "${RED}Outputs are different for ${name} on U = ${U}!${NC}\n"
			exit 1
		fi

	    echo -e "${name},${U},${t1},${t2},${t3},${t4}" >> results.csv
    done
done

printf "\n${GREEN}Finished running the experiment!\n${BLUE}Wrote algorithm times to results.csv${NC}\n"
