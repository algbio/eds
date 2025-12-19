#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

wget "https://www.cs.helsinki.fi/group/gsa/efg-mems/covid19-ecoli-efg.zip" --output-document=input/covid19-ecoli-efg.zip
unzip -p input/covid19-ecoli-efg.zip \
	inputs/covid19-100.fa \
	> input/covid19-100.fa

# change all ambiguous nucleotides to Ns
awk '{
	if (substr($0, 1, 1) == ">")
		{print}
	else
		{gsub(/[RYSWKMBDHV]/, "N"); print}
	}' input/covid19-100.fa > input/covid19-100-N.fa
