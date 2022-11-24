#!/bin/bash

for i in $(ls *corecov_coord.txt); do
    echo $i
    SAMPLES=$( perl filt_vcf_samples.pl $i freebayes_20220112_core.vcf )
    FILE_PREFIX="${i%_corecov_coord.txt}"
    OUTFILE=$FILE_PREFIX".vcf"
    /home/kirsten/Software/vcflib/bin/vcfkeepsamples temp.vcf $SAMPLES > $OUTFILE
done
    
