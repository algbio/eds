#!/bin/bash
set -e
set -o pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

viralmsa=$thisfolder/../ext/ViralMSA/ViralMSA.py
minimap2_folder=$thisfolder/../ext/minimap2
export PATH=$minimap2_folder:$PATH
datasets=datasets

for p in $viralmsa $minimap2_folder/minimap2 $datasets openssl
do
	if ! command -v $p >/dev/null 2>&1
	then
		echo "ERROR: $p not found!"; exit 1
	fi
done

# https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html#Random-sources
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

cd input

# 1. get reference
datasets download virus genome accession "NC_045512.2"
unzip -p ncbi_dataset.zip ncbi_dataset/data/genomic.fna > reference.fa
rm ncbi_dataset.zip

# 2. get datasets
for s in 100 1000 10000 100000
do
	gunzip -c SARS-CoV-2_2026-04-15.acc.gz | \
	       	shuf -n $s --random-source=<(get_seeded_random "eds") \
		> sample_$s.acc
	datasets download virus genome accession --inputfile sample_$s.acc
	unzip -p ncbi_dataset.zip ncbi_dataset/data/genomic.fna > covid_${s}_unaligned.fa
	rm ncbi_dataset.zip
	rm -Rf viralmsa_$s
	$viralmsa --sequences covid_${s}_unaligned.fa --reference reference.fa --output viralmsa_$s
	rm covid_${s}_unaligned.fa
	mv viralmsa_$s/covid_${s}_unaligned.fa.aln covid_${s}.fa
	bgzip --force covid_${s}.fa
done
