#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

mincard=$thisfolder/../../msa2eds-mincard
inputmsa=$thisfolder/input/msa.fa

mkdir output
cd output
ln -s $inputmsa msa.fa

for U in 4 8 16
do
	$mincard msa.fa $U
	mv msa.fa.eds.txt eds_U$U.txt
done
