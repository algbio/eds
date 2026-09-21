#!/bin/bash
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

seqtk=$([ -f "$thisfolder/../ext/seqtk/seqtk" ] && echo "$thisfolder/../ext/seqtk/seqtk" || echo "seqtk")

for p in wget xz $seqtk
do
	if ! command -v $p >/dev/null 2>&1
	then
		echo "ERROR: $p not found!"; exit 1
	fi
done

wget "https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Forschung/Projekte/seqana/MSA/chr19_100.aligned.fa.xz.00" --output-document=input/chr19_100.aligned.fa.xz.00
wget "https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Forschung/Projekte/seqana/MSA/chr19_100.aligned.fa.xz.01" --output-document=input/chr19_100.aligned.fa.xz.01
wget "https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Forschung/Projekte/seqana/MSA/chr19_100.aligned.fa.xz.02" --output-document=input/chr19_100.aligned.fa.xz.02

cat input/chr19_100.aligned.fa.xz.{00,01,02} | xz -d | $seqtk seq -U > input/chr19_100.aligned.uppercase.fa
