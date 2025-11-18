#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

datasets=datasets
iqtree=iqtree3

$datasets download genome accession GCF_000005845.2 \
	--filename input/e_coli_K-12_MG1655_dataset.zip
unzip -p input/e_coli_K-12_MG1655_dataset.zip \
	ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna \
	> input/e_coli_K-12_MG1655.fasta
$iqtree --alisim input/msa \
	--root-seq input/e_coli_K-12_MG1655.fasta,NC_000913.3 \
	--indel 0.01,0.01 -t input/tree.nwk \
	-af fasta \
	--seed 302288
