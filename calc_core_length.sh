#!/bin/bash

#The first awk command makes sure the start positions are inferior to the end positions. Then sort reorders the bed file using the start position.
#The second awk command merge the overlapping intervals.
#The third awk command computes the length of the intervals.

printf "Strain_id\tnew_id\tTot_length(bp)\n"

for f in *_core_filt.bed ; do
	str=$(awk -F '\t' 'NR==1 {print $1}' $f) ;
	awk -F '\t' '{
		if ($2<$3) {
			print $0
		} else {
			printf "%s\t%s\t%s\t%s\n",$1,$3,$2,$4
		}
	}' $f | sort -nk2 | awk -F '\t' 'NR==1 {
		st=$2;
		nd=$3
	}
	NR>1 {
		if ($2<nd+1) {
			nd=$3
		} else {
			print st "\t" nd;
			st=$2;
			nd=$3
		}
	}' | awk -F '\t' -v sdp="${f%_core_filt.bed}" -v strain=$str 'BEGIN {
			sum=0
		} 
		{
			sum=sum+$2-$1
		}
		END {
			print sdp "\t" strain "\t" sum
		}'
done
