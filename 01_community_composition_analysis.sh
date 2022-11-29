#!/bin/bash

#Description: This script will map the filtered reads against the core genes of the beebiome bacterial database, filter the mapping, use the mapping data to infer terminus coverage of each of the bacterial SDPs, use the qPCR data to infer bacterial loads of the samples, then use a R script to generate the figures. WARNING: this script will erase any existing mapping_full_db and community_composition directory and create new ones. NOTE: this script is voluntarily unoptimized to run more easily. You can easily make it faster by parallelizing the mapping and filtering of mapping data operations. Mos of the scripts here were written by Kirsten Ellegaard. Note that if you are using a subset of the samples or a subset of the database (or adding samples or strain genomes to the database), you will need to edit the perl "core_cov.pl" script in which the list of samples and genomes is hardcoded.

#Requirements:
#	* Having run the 00_data_preparation.sh script.
#	* Having downloaded the beebiome reference database and put it in a directory called "beebiome_db".
#	* bash, awk, perl, bwa, samtools, java, picard tools and R (with packages segmented, plyr, ape, data.table, ggplot2, gridExtra, RColorBrewer, reshape, ggrepel, qvalue and vegan) installed and in the path.

#Usage: ./01_community_composition_analysis.sh

if [ -d mapping_full_db ] ; then rm -r mapping_full_db ; fi
mkdir mapping_full_db #Create a QC directory if it does not exist yet.
if [ -d community_composition ] ; then rm -r community_composition ; fi
mkdir community_composition #Create a QC directory if it does not exist yet.

bwa index beebiome_db/beebiome_db

map_and_filter () {
	tmp="${1%_R1_*}"
    sample_tag="${tmp##*/}"
	file2="${tmp}_R2_combined.fastq.gz"

	echo "${sample_tag}"

	bwa mem -t 10 beebiome_db/beebiome_db $1 $file2 > "mapping_full_db/${sample_tag}_vs_db.sam"
	perl filter_sam_aln_length.pl "mapping_full_db/${sample_tag}_vs_db.sam" > "mapping_full_db/${sample_tag}_vs_db_filtered.sam"
	samtools view -bSu "mapping_full_db/${sample_tag}_vs_db_filtered.sam" | samtools sort -o "mapping_full_db/${sample_tag}_vs_db_filtered_sorted.bam"
	samtools view -F4 -h -@10 "mapping_full_db/${sample_tag}_vs_db_filtered_sorted.bam" > "mapping_full_db/${sample_tag}_tmp.sam"
	samtools view -F 0x800 -h "mapping_full_db/${sample_tag}_tmp.sam" | grep -E "NM:i:[0-4][[:blank:]]|^\@" > "mapping_full_db/${sample_tag}_edit_dist.sam"
	samtools view -bh "mapping_full_db/${sample_tag}_edit_dist.sam" | samtools sort -o "mapping_full_db/${sample_tag}_vs_db_filt2x_sorted.bam"
	samtools index "mapping_full_db/${sample_tag}_vs_db_filt2x_sorted.bam"
	rm "mapping_full_db/${sample_tag}_vs_db_filtered.sam"
	rm "mapping_full_db/${sample_tag}_vs_db_filtered_sorted.bam"
	rm "mapping_full_db/${sample_tag}_tmp.sam"
	rm "mapping_full_db/${sample_tag}_edit_dist.sam"
	
	perl filter_sam_aln_length_low "mapping_full_db/${sample_tag}_vs_db.sam" > "mapping_full_db/${sample_tag}_vs_db_unmapped.sam"
	samtools view -bSu "mapping_full_db/${sample_tag}_vs_db_unmapped.sam" | samtools sort -o "mapping_full_db/${sample_tag}_vs_db_unmapped_sorted.bam"
	rm "mapping_full_db/${sample_tag}_vs_db_unmapped.sam"
	samtools view -bh "mapping_full_db/${sample_tag}_vs_db.sam" > "mapping_full_db/${sample_tag}_vs_db.bam"
	rm "mapping_full_db/${sample_tag}_vs_db.sam"
	java -Xmx2g -jar ~/Softwares/picard.jar SamToFastq I="mapping_full_db/${sample_tag}_vs_db_unmapped_sorted.bam" FASTQ="mapping_full_db/${sample_tag}_vs_db_unmapped_R1.fastq" SECOND_END_FASTQ="mapping_full_db/${sample_tag}_vs_db_unmapped_R2.fastq" VALIDATION_STRINGENCY=SILENT #Picard tool converts the bam files to fastq files (forward and reverse reads).
	gzip "mapping_full_db/${sample_tag}_vs_db_unmapped_R1.fastq"
	gzip "mapping_full_db/${sample_tag}_vs_db_unmapped_R2.fastq"
	rm "mapping_full_db/${sample_tag}_vs_db_unmapped_sorted.bam"
}

terminus_coverage_regression () {
	phylotype="${1%_single_ortho.txt}"
	perl core_cov.pl $1
	Rscript core_cov.R "${phylotype}_corecov.txt"
}

echo "Starting mapping and filtering."

file_list=($(ls filtered_reads/*_R1_combined.fastq.gz))

for f in "${file_list[@]}" ; do 
	map_and_filter $f
done

echo "Mapping and filtering done. Starting terminus coverage regression."

cd community_composition
ln -s ../core_cov.pl .
ln -s ../core_cov.R .
ln -s ../mapping_full_db/*_filt2x_sorted.bam .
ln -s ../beebiome_db/bed_files/*.bed .
ln -s ../beebiome_db/single_ortho_files/*_single_ortho.txt .

phylo_list=($(ls *_single_ortho.txt))

for f in "${phylo_list[@]}" ; do
	terminus_coverage_regression $f
	unlink $f
done

for f in *.bam ; do unlink $f ; done
for f in *.bed ; do unlink $f ; done

cd ..

echo "Terminus overage regression done. Starting to produce the community composition figures."

Rscript figures_community_composition.R

echo "Community composition figures done."

