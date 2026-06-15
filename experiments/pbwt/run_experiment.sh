#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd "$thisfolder"

mincard_naive="$thisfolder/mincard-naive"
mincard_recursive="$thisfolder/mincard-recursive"
mincard_rmq="$thisfolder/mincard-rmq"
inputmsa="$thisfolder/input/covid_100000.fa.gz"
usrbintimeformat="%e"

outdir="output"
mkdir -p "$outdir"
rm -f "$outdir"/*
cd "$outdir"
ln -s "$inputmsa" msa.fa

U_values=(4 8 16 32 64 128)
printf "Running mincard with pBWT and different range max strategies (naive, recursive, rmq) on different U values %s \n" "${U_values[*]}"

tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

declare -a naive_times
declare -a recursive_times
declare -a rmq_times

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m'

# mincard with pbwt (gaps as symbols strategy)
for U in "${U_values[@]}"
do
    printf "\n${BLUE}Comparison for U = $U${NC}\n"
	# naive
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard_naive" msa.fa --gaps-as-symbols --pbwt -U $U -o mincard_naive_U${U}.eds
	t1=$(<"$tmpfile")
	naive_times+=("$t1")
	# recursive
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard_recursive" msa.fa --gaps-as-symbols --pbwt -U $U -o mincard_recursive_U${U}.eds
	t2=$(<"$tmpfile")
	recursive_times+=("$t2")
	# rmq
	/usr/bin/time -f"$usrbintimeformat" -o "$tmpfile" "$mincard_rmq" msa.fa --gaps-as-symbols --pbwt -U $U -o mincard_rmq_U${U}.eds
	t3=$(<"$tmpfile")
	rmq_times+=("$t3")

	if cmp -s mincard_naive_U${U}.eds mincard_recursive_U${U}.eds &&
       cmp -s mincard_naive_U${U}.eds mincard_rmq_U${U}.eds;
	then
		printf "${GREEN}Outputs are identical for U = ${U}!${NC}\n"
		printf "${PURPLE}pBWT max range: naive time %.3fs vs. recursive time %.3fs vs. rmq time %.3fs${NC}\n" "$t1" "$t2" "$t3"
	else 
		printf "${RED}Outputs are different for U = ${U}!${NC}\n"
		exit 1
	fi
done
# Write results to table
echo -e "U value,Naive,Recursive,RMQ" > results.csv
for i in "${!U_values[@]}";
do
	echo -e "${U_values[$i]},${naive_times[$i]},${recursive_times[$i]},${rmq_times[$i]}" >> results.csv
done
printf "\n${BLUE}Wrote algorithm times to results.csv${NC}\n"