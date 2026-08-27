#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$thisfolder"
export LC_NUMERIC="en_US.UTF-8"

mincard="$thisfolder/../../mincard"
minsize="$thisfolder/../../minsize"

# MSA fasta files from the 3 datasets
input_ecoli_msa="$thisfolder/../e_coli_sim/input/msa.fa"
input_covid_msa="$thisfolder/../covid/input/covid_100000.fa.gz"
input_chr19_msa="$thisfolder/../chr19/input/chr19_100.aligned.uppercase.fa"

# Create links to msa files
outdir="output"
mkdir -p "$outdir"
rm -f "$outdir"/*
cd "$outdir"

ln -s "$input_ecoli_msa" ecoli_msa.fa
ln -s "$input_covid_msa" covid_msa.fa
ln -s "$input_chr19_msa" chr19_msa.fa

# Keep them in an array for easy iteration
declare -A msa

datasets=(
    ecoli
    covid
    chr19
)

msa[ecoli]=ecoli_msa.fa
msa[covid]=covid_msa.fa
msa[chr19]=chr19_msa.fa

# Store algorithm flags
declare -A flags
declare -A alg

algos=(
	mincard_gaps
	mincard_gapless
	minsize_gaps
	minsize_gapless
)

flags[mincard_gaps]="--gaps"
flags[minsize_gaps]="--gaps"
flags[mincard_gapless]=""
flags[minsize_gapless]=""

alg[mincard_gaps]="$mincard"
alg[mincard_gapless]="$mincard"
alg[minsize_gaps]="$minsize"
alg[minsize_gapless]="$minsize"

# Create temporary files for collecting results
outfile=$(mktemp)
timefile=$(mktemp)
trap 'rm -f "$outfile"' EXIT
trap 'rm -f "$timefile"' EXIT

usrbintimeformat="%e"

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m'

# Iterate through different [L, U] bounds

values=("1 16" "16 32" "32 64")
printf "Running mincard/minsize with gaps/gapless strategies on MSAs with LU bounds"
for pair in "${values[@]}"; do
	printf " [%s]" "${pair/ /, }"
done
printf "\n"


# Write results to csv table
echo -e "Dataset,Algorithm,L,U,Cardinality,Size,Avg. Segment,Time" > results.csv

for dataset in "${datasets[@]}"; do
	input_msa="${msa[$dataset]}"
	
	for pair in "${values[@]}"; do
		read lower upper <<< "$pair"
		printf "\nDataset: ${BLUE}${dataset}${NC} |  L = ${BLUE}${lower}${NC} U = ${BLUE}${upper}${NC}\n"
		
		for algo in "${algos[@]}"; do
		  /usr/bin/time -f"$usrbintimeformat" -o "$timefile" \
				${alg[$algo]} "$input_msa" \
   		  -L "$lower" -U "$upper" \
  		  -q --stats \
   		  ${flags[$algo]} > "$outfile" 2>&1;
			
			time=$(cat "$timefile")
			read card size <<< "$(head -n1 "$outfile")"
			read min max avg <<< "$(tail -n1 "$outfile")"

			printf "Algorithm: ${GREEN}${algo}${NC} | cardinality = ${PURPLE}${card}${NC} | gap-aware size = ${PURPLE}${size}${NC} | segment sizes: min = ${PURPLE}${min}${NC}, max = ${PURPLE}${max}${NC}, avg = ${PURPLE}${avg}${NC} | ${BLUE}${time}s${NC}\n"
	    echo -e "${dataset},${algo},${lower},${upper},${card},${size},${avg},${time}" >> results.csv
		done
	done
done

printf "\n${GREEN}Finished running the experiment!\n"
printf "${BLUE}Wrote algorithm times and stats to results.csv${NC}\n"
