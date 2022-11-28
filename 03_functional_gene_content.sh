#! /bin/bash

#Description: This script will assemble the filtered reads which did not map the host genome (so it first needs to map the unmapped read from step 01 to the host genome, and combine the new unmapped reads with the reads mapping the beebiome_db), filter them, detect the open reading frames, filter them, map the reads against the orfs database, filter the mappings and extract the coverage. On another hand, the orfs are blasted against the beebiome_db to bin the contigs to the SDPs. The orfs are also assigned to OGs per SDP. Finally, a R script generates all the figures. WARNING: this script will erase any existing assemblies and functional_gene_content directory and create new ones. NOTE: this script is not optimized and will run slowly. You can easily make it faster. Most of the scripts were written by Kirsten Ellegaard.

#Requirements:
#	* Having run successfully the 00_data_preparation.sh script.
#	* Having downloaded the beebiome reference database and put it in a directory called "beebiome_db".
#	* bash, awk, perl, python, SPAdes, bwa, samtools, java, picard tools, prodigal, othofinder, eggnog and R (with packages ape, data.table, ggplot2, gridExtra, RColorBrewer, reshape, ggrepel, qvalue and scales) installed and in the path.

#Usage ./03_functional_gene_content.sh

if [ -d assemblies ] ; then rm -r assemblies ; fi
mkdir assemblies
if [ -d functional_gene_content ] ; then rm -r functional_gene_content ; fi
mkdir functional_gene_content

#First, map beebiome-unmapped reads to the honeybee genome, and extract honeybee-unmapped reads.

bwa index beebiome_db/honeybee_genome

map_and_filter () {
	tmp="${1%_vs_db_unmapped_R1*}"
    	sample_tag="${tmp##*/}"
	file2="${tmp}_vs_db_unmapped_R2.fastq.gz"

	echo "${sample_tag}"

	bwa mem -t 10 beebiome_db/honeybee_genome.fasta $1 $file2 > "mapping_full_db/${sample_tag}_unmapped_vs_host.sam"
	perl filter_sam_aln_length_low.pl "mapping_full_db/${sample_tag}_unmapped_vs_host.sam" > "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped.sam"
	samtools view -bSu "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped.sam" | samtools sort -o "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_sorted.bam"
	rm "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped.sam"
	samtools view -bh "mapping_full_db/${sample_tag}_unmapped_vs_host.sam" > "mapping_full_db/${sample_tag}_unmapped_vs_host.bam"
	rm "mapping_full_db/${sample_tag}_unmapped_vs_host.sam"
	java -Xmx2g -jar ~/Softwares/picard.jar SamToFastq I="mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_sorted.bam" FASTQ="mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_R1.fastq" SECOND_END_FASTQ="mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_R2.fastq" VALIDATION_STRINGENCY=SILENT #Picard tool converts the bam files to fastq files (forward and reverse reads).
	gzip "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_R1.fastq"
	gzip "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_R2.fastq"
	rm "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_sorted.bam"
	cat "filtered_reads/${sample_tag}_R1_combined.fastq.gz" "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_R1.fastq.gz" > "assemblies/${sample_tag}_non-host_R1.fastq.gz"
	cat "filtered_reads/${sample_tag}_R2_combined.fastq.gz" "mapping_full_db/${sample_tag}_unmapped_vs_host_unmapped_R2.fastq.gz" > "assemblies/${sample_tag}_non-host_R2.fastq.gz"
}

assemble () {
	tmp="${1%_non-host_R1*}"
        sample_tag="${tmp##*/}"
        file2="${tmp}_non-host_R2.fastq.gz"

	gunzip $1
	gunzip $file2

	python /home/kirsten/Software/SPAdes-3.10.1-Linux/bin/spades.py --meta -1 "${1%.gz}"  -2 "${file2%.gz}" -t 10 -o "assemblies/${sample_tag}"
	mv "assemblies/${sample_tag}/contigs.fasta" "assemblies/${sample_tag}_contigs.fasta"
	rm -r "assemblies/${sample_tag}"

	gzip $1
	gzip $file2
}

filter_contigs_predict_orfs_filter_orfs () {
	awk -F'_' 'BEGIN { k=0 } $1~">NODE" { if ($4 >= 300 && $6 >= 1) { k=1 } else { k=0 } } { if (k==1) { print $0 } }' $1 > "${1%.fasta}_filtered.fasta"
	prodigal -i "${1%.fasta}_filtered.fasta" -d "${1%_contigs*}.ffn" -a "${1%_contigs*}.faa" -m -p meta
	perl filt_ffn.pl "${1%_contigs*}.ffn"
}

map_against_orfs_filter_and_extract_coverage () {
	tmp="${1%_non-host_R1*}"
        sample_tag="${tmp##*/}"
        file2="${tmp}_non-host_R2.fastq.gz"
	
	echo "${sample_tag}"

        bwa mem -t 10 "functional_gene_content/orfs_and_genes_catalogue.ffn" $1 $file2 > "functional_gene_content/${sample_tag}_non-host_vs_catalogue.sam"
	perl filter_sam_aln_length.pl "functional_gene_content/${sample_tag}_non-host_vs_catalogue.sam" > "functional_gene_content/${sample_tag}_non-host_vs_catalogue_filtered.sam"
        samtools view -bSu "functional_gene_content/${sample_tag}_non-host_vs_catalogue_filtered.sam" | samtools view -bh -F 0x800 | samtools sort -o "functional_gene_content/${sample_tag}_non-host_vs_catalogue_filtered_sorted.bam"

	samtools coverage "functional_gene_content/${sample_tag}_non-host_vs_catalogue_filtered_sorted.bam" > "functional_gene_content/${sample_tag}_catalogue_cov.txt"
} 

bwa index beebiome_db/honeybee_genome.fasta

echo "Preparing the reads."

file_list=($(ls mapping_full_db/*_vs_db_unmapped_R1.fastq.gz))

for f in "${file_list[@]}" ; do map_and_filter $f ; done

echo "Reads prepared. Starting assemblies."

file_list=($(ls assemblies/*_non-host_R1.fastq.gz))

for f in "${file_list[@]}" ; do assemble $f ; done

echo "Assemblies done. Starting filtering, orf prediction and orf filtering."

file_list=($(ls assemblies/*_contigs.fasta))

for f in "${file_list[@]}" ; do filter_contigs_predict_orfs_filter_orfs $f ; done

echo "ORF predicted and filtered. Starting binning."

cat assemblies/*_filt.ffn > functional_gene_content/all_filt.ffn
cat assemblies/*_filt.faa > functional_gene_content/all_filt.faa

makeblastdb -dbtype 'nucl' -in beebiome_db/all_beebiome_genes_concat.ffn
bwa index beebiome_db/all_beebiome_genes_concat.ffn

blastn -db beebiome_db/all_beebiome_genes_concat.ffn -query functional_gene_content/all_filt.ffn  -num_threads 10 -max_target_seqs 1 -outfmt '6 qseqid sseqid evalue bitscore pident length qcovs' > functional_gene_content/all_filt_ffn.blastn

awk -F'\t' ' $5>=50 && $7>=50 && $3<1e-05 { pid[$1][$2]=$5; lin[$1][$2]=$0 } END { for (i in pid) { bestpid=0; best=0; for (j in pid[i]) { if (pid[i][j]>bestpid) { bestpid=pid[i][j]; best=j } } print lin[i][best] } }' functional_gene_content/all_filt_ffn.blastn > functional_gene_content/all_filt_ffn.blastn.1.filt #Filtering blast results

awk -F'\t' ' NR==FNR { a[$1]=$4 ; b[$1]=$3 } NR!=FNR { split($1,orf,/_/) ; split($2,x,/_/) ; y=x[1] ; printf "%s\t%s_%s_%s_%s_%s_%s_%s\t%s\t%s\t%s\t%s\n",$1,orf[1],orf[2],orf[3],orf[4],orf[5],orf[6],orf[7],$2,x[1],b[y],a[y] } ' db_metadata_GB.txt functional_gene_content/all_filt_ffn.blastn.1.filt > functional_gene_content/all_filt_ffn.blastn.1.filt.SDP #Adding SDP information.

grep '^>' functional_gene_content/all_filt.ffn | sed 's/>//g' > all_filt_ORFs_list.txt

awk -F'_' ' { printf "%s\t%s_%s_%s_%s_%s_%s_%s\n",$0,$1,$2,$3,$4,$5,$6,$7 } ' all_filt_ORFs_list.txt | awk -F'\t' ' NR==FNR { c[$2]++ } NR!=FNR { a[$2,$6]++ ; b[$2]=$2 ; d[$6]=$6 } END { printf "contig" ; for (l in d) { printf "\t%s",d[l] } for (k in b) { printf "\n%s",b[k] ; for ( l in d) { if ( a[b[k],d[l]]=="" ) { a[b[k],d[l]]=0 } printf "\t%.4f",a[b[k],d[l]]/c[k] } } } ' - functional_gene_content/all_filt_ffn.blastn.1.filt.SDP > functional_gene_content/contig_x_SDP.txt #Table of proportions of orfs belonging to each SDP for all contigs.

awk -F'\t' ' NR==1 { for (i=2;i<=NF;i++) { sdp[i]=$i } } NR>1 && NR==FNR { imax=0 ; max=0 ; cum=0 ; for (i=2;i<=NF;i++) { cum=cum+$i ; if ( $i>max ) { imax=i ; max=$i } } ; if ( max>=(0.8*cum) && cum>=0.1 ) { contigl[$1]=$1 ; contig[$1]=(sdp[imax] FS max FS cum) } } NR>FNR { split($1,a,"_") ; x=a[1]"_"a[2]"_"a[3]"_"a[4]"_"a[5]"_"a[6]"_"a[7] ; if (x in contigl) { print $1 FS x FS contig[x] } else { print $1 FS x FS "unassigned" FS "0" FS "0" } } ' functional_gene_content/contig_x_SDP.txt functional_gene_content/all_filt_ORFs_list.txt > functional_gene_content/all_filt_ORFs_list_assigned.txt

echo "Binning done. Starting to map the reads back to the gene and orf catalogue."

#Running orthofinder to define orthogroups per SDP.
mkdir functional_gene_content/faa_files
for f in beebiome_db/faa_files/*.faa ; do
	tmp=${f%.faa}
	loc_tag=${tmp##*/}
	sdp=$(awk -F'\t' -v file=${loc_tag} ' $1~file {print $4} ' db_metadata_GB.txt)
	mkdir -p "./functional_gene_content/faa_files/${sdp}"
	cat $f > "./functional_gene_content/faa_files/${sdp}/${sdp}_${loc_tag}"
done
mkdir -p "unassigned"
for f in assemblies/*_filt.faa ; do
	tmp=${f%_filt.faa}
	sample="${tmp##*/}"
	awk -v smpl=$sample -F'\t' 'NR==FNR { a[$1]=$6 } NR>FNR && $1 ~ /^>/ { orf=substr($1,2) } NR>FNR { print $0 >> "./functional_gene_content/faa_files/"a[orf]"/"a[orf]"_"smpl }' ../all_filt_ORFs_list_assigned.txt $f
done
> functional_gene_content/all_filt_ffn.blastn.1.filt.SDP.OG
for d in functional_gene_content/faa_files/* ; do
	sdp=${d##*/}
	echo $sdp
	python orthofinder.py -og -t 10 -f $d
	tmp=$(ls $d/OrthoFinder/)
	sed 's/://g' $d/OrthoFinder/$tmp/Orthogroups/Orthogroups.txt | awk -v sdp=$sdp ' NR==FNR { a[$1]=$0 } NR>FNR { for (i=2;i<=NF;i++) { b[$i]="\t"sdp"_"$1 } } END { for (x in a) { printf "%s%s\n",a[x],b[x] } }' functional_gene_content/all_filt_ffn.blastn.1.filt.SDP - | awk -F'\t' 'NF>6' - >> functional_gene_content/all_filt_ffn.blastn.1.filt.SDP.OG
	sed 's/://g' $d/OrthoFinder/$tmp/Orthogroups/Orthogroups.txt | awk -v sdp=$sdp ' { print sdp"_"$0 } ' - >> functional_gene_content/all_OGs.txt
done
> functional_gene_content/core_OGs_list.txt
for f in beebiome_db/*_ortho_filt.txt ; do
        tmp="${f%_single_ortho_filt.txt}"
	phylo=${f##*/}
        echo $phylo
        awk -v x=$phylo ' NR==FNR { for (i=2;i<=NF;i++) { a[$i]=$1 } } NR>FNR { for (i=2;i<=NF;i++) { printf "%s\t%s\t%s\n",a[$i],$i,x } } ' functional_gene_content/all_OGs.txt $f >> functional_gene_content/core_OGs_list.txt
done #Looking for coreness information for each OG.

#Creating the genes and orfs catalogue for the mapping.

cat functional_gene_content/all_filt.ffn beebiome_db/all_beebiome_genes_concat.ffn > functional_gene_content/orfs_and_genes_catalogue.ffn
bwa index functional_gene_content/orfs_and_genes_catalogue.ffn

for f in functional_gene_content/*_non-host_R1.fastq.gz ; do 
	map_against_orfs_filter_and_extract_coverage $f
done

#Using eggnog to annotate the database.
cat beebiome_db/faa_files/*.faa functional_gene_content/all_filt.faa > functional_gene_content/orfs_and_genes_catalogue.faa
python emapper.py -m diamond --no_annot --no_file_comments --cpu 10 -i functional_gene_content/orfs_and_genes_catalogue.faa -o functional_gene_content/all_ogs_eggnog
python emapper.py --annotate_hits_table functional_gene_content/all_ogs_eggnog.emapper.seed_orthologs --no_file_comments -o signifORFs_eggnog_annot --cpu 10

#Creating the figures
echo "Generating the figures."
Rscript figures_functional_gene_content.R

echo "Functional gene content analysis done."

