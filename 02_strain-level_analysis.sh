#!/bin/bash

#Description: This script will map the filtered reads against the core genes of the beebiome bacterial database reduced to one genome per SDP. It will then filter the alignments, use the mapping data with freebayes to find the Single nucleotide Variants in each sample. From this, it creates statistic files of the variability, but also cumulative curves of the total cumulative proportion of variable sites. Finally, a R script generates all the figures. WARNING: this script will erase any existing mapping_red_db and strain-level_analysis directory and create new ones. NOTE: this script is not optimized and will run slowly. You can easily make it faster. Most of the scripts were written by Kirsten Ellegaard. 

#Requirements:
#	* Having run successfully the 00_data_preparation.sh script.
#	* Having downloaded the beebiome reference database and put it in a directory called "beebiome_db".
#	* bash, awk, perl, bwa, samtools, java, picard tools, vcflib, freebayes and R (with packages ape, data.table, ggplot2, gridExtra, RColorBrewer, reshape, ggrepel, qvalue, scales and PERMANOVA) installed and in the path.

#Usage ./02_strain-level_analysis.sh

if [ -d mapping_red_db ] ; then rm -r mapping_red_db ; fi
mkdir mapping_red_db 
if [ -d strain-level_analysis ] ; then rm -r strain-level_analysis ; fi
mkdir strain-level_analysis

bwa index beebiome_db/beebiome_red_db

map_and_filter () {
        tmp="${1%_R1_*}"
        sample_tag="${tmp##*/}"
        file2="${tmp}_R2_combined.fastq.gz"

        echo "${sample_tag}"

        bwa mem -t 10 beebiome_db/beebiome_red_db $1 $file2 > "mapping_red_db/${sample_tag}_vs_db.sam"
        perl filter_sam_aln_length.pl "mapping_red_db/${sample_tag}_vs_db.sam" > "mapping_red_db/${sample_tag}_vs_db_filtered.sam"
        samtools view -bSu "mapping_red_db/${sample_tag}_vs_db_filtered.sam" | samtools sort -o "mapping_red_db/${sample_tag}_vs_db_filtered_sorted.bam"
        samtools view -F4 -h -@10 "mapping_red_db/${sample_tag}_vs_db_filtered_sorted.bam" > "mapping_red_db/${sample_tag}_tmp.sam"
        samtools view -F 0x800 -h "mapping_red_db/${sample_tag}_tmp.sam" | grep -E "NM:i:[0-4][[:blank:]]|^\@" > "mapping_red_db/${sample_tag}_edit_dist.sam"
        samtools view -bh "mapping_red_db/${sample_tag}_edit_dist.sam" | samtools sort -o "mapping_red_db/${sample_tag}_vs_db_filt2x_sorted.bam"
        rm "mapping_red_db/${sample_tag}_vs_db_filtered.sam"
        rm "mapping_red_db/${sample_tag}_vs_db_filtered_sorted.bam"
        rm "mapping_red_db/${sample_tag}_tmp.sam"
        rm "mapping_red_db/${sample_tag}_edit_dist.sam"
        samtools view -bh "mapping_red_db/${sample_tag}_vs_db.sam" > "mapping_red_db/${sample_tag}_vs_db.bam"
        rm "mapping_red_db/${sample_tag}_vs_db.sam"
	java -jar picard.jar AddOrReplaceReadGroups I="mapping_red_db/${sample_tag}_vs_db_filt2x_sorted.bam" O="mapping_red_db/${sample_tag}_rg.bam" RGID=${sample_tag} RGLB=lib1 RGPL=illumina RGPU=none RGSM=${sample_tag}
	samtools index "mapping_red_db/${sample_tag}_rg.bam"
	rm "mapping_red_db/${sample_tag}_vs_db_filt2x_sorted.bam"
}

snv_calling () {
	
}

echo "Starting mapping and filtering."

file_list=($(ls filtered_reads/*_R1_combined.fastq.gz))

for f in "${file_list[@]}" ; do
        map_and_filter $f
done

echo "Mapping and filtering done. Starting SNV calling."

cd strain-level_analysis
ln -s ../core_cov_red.pl .
ln -s ../core_cov.R .
ln -s ../mapping_red_db/*_rg.bam .
ln -s ../beebiome_db/red_bed_files/*.bed .
ln -s ../beebiome_db/single_ortho_files/*_single_ortho.txt .

for f in *_rg.bam ; do echo $f >> bamfile_list.txt ; done 
perl core_cov_red.pl --bam-list bamfile_list.txt
for f in *_corecov.txt ; do Rscript core_cov.R $f ; done

ln -s ../filt_core_bed.pl .
for i in *coord.txt ; do perl filt_core_bed.pl $i ; done
cat *filt.bed > all_filt.bed

#ln -s ../calc_filt_core_length.pl .
#perl calc_filt_core_length.pl > table_core_length.txt
./calc_core_length.sh > table_core_length.txt

freebayes-parallel <(fasta_generate_regions.py beebiome_red_db.fai 100000) 30 -f beebiome_red_db -C 5 --min-alternate-fraction 0.1 --pooled-continuous --min-coverage 10 -i -X -u --bam-list bamfile_list.txt >  freebayes.vcf

vcfintersect --bed all_filt.bed freebayes.vcf | vcfbreakmulti > freebayes_core.vcf

./filter_vcf_samples.sh

rm freebayes_core.vcf
for i in *.vcf ; do echo $i ; perl filter_snvs.pl $i >> summary_filtering.txt ; done

for i in *filt.freq ; do echo $i ; SDP=${i%_filt.freq};COORD=$SDP"_corecov_coord.txt"; perl summarize_snps_host.pl table_core_length.txt $i $COORD; done

cat *tot_var.txt >> tot_var_all.txt
cat *sample_var_host.txt > all_sample_var_host.txt

for i in *filt.freq ; do echo $i ; perl cum_curve_host.pl table_core_length.txt $i ; done

for i in *filt.freq ; do echo $i; perl calc_jaccard.pl $i; done
for i in *shared_fraction.txt ; do PREFIX=${i%_filt_shared*} ; OUTFILE=$PREFIX"_dist_matrix.txt" ; Rscript distance_matrix.R $i $OUTFILE ; done

Rscript figures_strain-level_analyses.R


