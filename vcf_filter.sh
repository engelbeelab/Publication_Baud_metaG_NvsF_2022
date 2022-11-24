#!/bin/bash

/home/kirsten/Software/vcflib/bin/vcfintersect --bed all_filt.bed freebayes_20201215_s.vcf | /home/kirsten/Software/vcflib/bin/vcfbreakmulti > freebayes_20201215_s_core.vcf
