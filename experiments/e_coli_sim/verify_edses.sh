#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder/output

junctions=$thisfolder/../ext/junctions/bin/junctions
threads=8; if [ $# -gt 0 ] ; then threads=$1 ; fi

# remove gaps and split MSA into single files
awk '{
	if (substr($0, 1, 1) == ">")
		{print}
	else
		{gsub(/-/, ""); print}
	}' msa.fa | \
	awk 'BEGIN {RS=">"; FS="\n"}
		NR>1 {
			$1="";
			n += 1;
			outfile = "sequence-" n "-nogaps.eds";
			print $0 > outfile;
			close(outfile);
		}'

# check for intersections
for eds in *.eds
do
	echo "Verifying $eds..."
	ls sequence-* | \
		xargs -P $threads -I % \
		$junctions intersect $eds %
	echo "done."
done
