#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

seqtk=seqtk

wget "https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Forschung/Projekte/seqana/MSA/chr19_100.aligned.fa.xz.00" --output-document=input/chr19_100.aligned.fa.xz.00
wget "https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Forschung/Projekte/seqana/MSA/chr19_100.aligned.fa.xz.01" --output-document=input/chr19_100.aligned.fa.xz.01
wget "https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Forschung/Projekte/seqana/MSA/chr19_100.aligned.fa.xz.02" --output-document=input/chr19_100.aligned.fa.xz.02

cat input/chr19_100.aligned.fa.xz.{00,01,02} | xz -d | $seqtk seq -U > input/chr19_100.aligned.uppercase.fa
