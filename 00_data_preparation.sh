#!/bin/bash

#Description: This script will trim, QC, and merge the raw_reads into the filtered reads. This step is necessary for all downstream operations. WARNING: this script will erase any existing QC and filtered_reads directory and create a new one. NOTE: this script is voluntarily unoptimized to run more easily. You can easily make it much faster by running fastqc in the background (you do not need it to be finished before running trimmomatic or the next iteration of the loop) and parallelize the loop runs. If you are using a cluster, consider running all trimmomatic instances in parallel (same for fastqc), and then do the merge (start an istance for each sample name).

#Requirements: 
#	* A directory called "raw_reads" with all the raw reads (Illumina PE150 prepared with Nextera adapters) with the names in the format: "${sample_name}_L${sequencing_lane_number}_R${orientation_1_or_2}_${unique_tag}.fastq.gz"
#	* The NexteraPE-PE.fa file, containing the Nextera PE adapter sequences.
#	* bash, awk, fastqc, Trimmomatic-0.35, java installed and in the path.

#Usage: ./00_data_preparation.sh

if [ -d QC ] ; then rm -r QC ; fi
mkdir QC #Create a QC directory if it does not exist yet.
if [ -d filtered_reads ] ; then rm -r filtered_reads ; fi
mkdir filtered_reads #Create a QC directory if it does not exist yet.

filter_and_QC () {
	tmp="${1%_R1_*}"
	sample_tag="${tmp##*/}"

	file2=$(find raw_reads/ -maxdepth 1 -name "${sample_tag}_R2_*")

	fastqc -o QC/ --noextract $1 $file2
	
	java -jar trimmomatic-0.35.jar PE -threads 10 "${1}" "${file2}" "filtered_reads/${sample_tag}_R1_paired.fastq.gz" "filtered_reads/${sample_tag}_R1_unpaired.fastq.gz" "filtered_reads/${sample_tag}_R2_paired.fastq.gz" "filtered_reads/${sample_tag}_R2_unpaired.fastq.gz" ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 #Trimming of the reads: removing the adapters, the short reads, etc.

	rm "filtered_reads/*_unpaired.fastq.gz"

	fastqc -o QC/ --noextract "filtered_reads/${sample_tag}_R1_paired.fastq.gz" "filtered_reads/${sample_tag}_R2_paired.fastq.gz"

}

merge_filt_reads () {
	tmp="${1%_R1_*}"
        sample_tag="${tmp##*/}"
        sample_name="${sample_tag%_L[0-9]*}"

	> "${tmp}_R1_combined.fastq.gz"
        > "${tmp}_R2_combined.fastq.gz"

	for f in filtered_reads/"${sample_name}"_*_R1_paired.fastq.gz ; do
		R2_file="${f2%_R1*}_R2_paired.fastq.gz"
        	cat $f >> "${tmp0}_R1_combined.fastq.gz"
        	cat ${R2_file} >> "${tmp0}_R2_combined.fastq.gz"
	done
	rm ${tmp}*
}

echo "Filtering and QC."

for f in raw_reads/*_R1_*.fasq.gz ; do
	echo $f
	filter_and_QC $f
done

echo "Merging filtered reads."

file_list=($(ls filtered_reads/*_paired.fastq.gz | awk -F '_' -v OFS='_' '{print $1 FS $2}' | awk '!a[$0]++'))

for f in "${file_list[@]}" ; do
	echo $f
	merge_filt_reads $f
done

echo "Done."
