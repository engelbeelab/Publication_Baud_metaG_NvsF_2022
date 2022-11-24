# bee_gut_metagen_NF

Scripts for the "Ecological turnover of strain-level diversity in the honeybee gut microbiota between nurses and foragers" manuscript.
=======

This repository contains the necessary scripts and metadata required to reproduce the analyses and figures from the "Ecological turnover of strain-level diversity in the honeybee gut microbiota between nurses and foragers" manuscript. They use metagenomic data from honeybee gut DNA sequencing with Illumina.

Most of the scripts and analyses from this manuscript have been adapted, developed or taken directly from the work of Kirsten Ellegaard, in particular the following publications:

> Kirsten Maren Ellegaard & Philipp Engel. **Genomic diversity landscape of the honey bee gut microbiota**; _Nature Communications_ **10**, Article number: 446 (2019).
> PMID: 30683856;
> doi:[10.1038/s41467-019-08303-0](https://www.nature.com/articles/s41467-019-08303-0)

> Kirsten Maren Ellegaard, Shota Suenami, Ryo Miyasaki, Philipp Engel. **Vast differences in strain-level diversity in the gut microbiota of two closely related honey bee species**; _Current Biology_ **10**, Epub 2020 Jun 11.
> PMID: 32531278;
> doi: [10.1016/j.cub.2020.04.070](https://www.cell.com/current-biology/fulltext/S0960-9822(20)30586-8)

If you are using the scripts found here, please cite them as well.


 
Pre-requisites
----------

The scripts and analyses rely on a genomic database built by Kirsten Ellegaard ([zenodo_link](https://zenodo.org/record/4661061#.YGmkRy0RoRA)), which contains mostly genomes of bacteria previously isolated from guts of *Apis mellifera*, *Apis cerana* and bumble bee species.

In terms of software, the scripts require:

* Python 3 (version 3.6 or higher)
* Bash
* R  (version 4.2.1, including packages: "segmented", "plyr", "ape" (version 5.6-2), "data.table" (version 1.14.2), "qvalue" (version 2.28.0), "PERMANOVA" (version 0.2.0), "vegan" (2.6-2), "ggplot2" (version 3.3.6), "reshape" (version 0.8.9), "gridExtra" (version 2.3), "ggrepel" (version 0.9.1), "RColorBrewer" (version 1.1-3) and "scales" (version 1.2.1))
* [samtools](http://www.htslib.org) (version  1.15-8-gbdc5bb8, should be added to the path)
* [bwa](https://github.com/lh3/bwa) (for the mapping) (version 0.7.15-r1142-dirty, should be added to the path)
* [SPAdes] (https://github.com/ablab/spades) (for the assemblies) (version 3.10.1, should be added to the path)
* [freebayes] (https://github.com/freebayes/freebayes) (for the detection of single nucleotide polymorphisms) (version v1.3.1-19-g54bf409, should be added to the path)
* [fastQC] (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (quality check of the input data) (should be added to the path)
* [Trimmomatic] (https://github.com/usadellab/Trimmomatic) (trimming of the reads and removal of the sequencing adapters) (version 0.35, should be added to the path)


Running the pipeline
--------

**Data preparation: mapping and filtering**

First, download all he necessary data. These are the raw sequencing reads (which should be put in a `raw_reads` directory), the reference database (`beebiome_db`) and the qPCR data (`qPCR_data`).
Then you can run `00_data_preparation.sh` script. \
```./00_data_preparation.sh```

**Community composition analysis**

This one should be straightforward. \

```./01_community_composition_analysis.sh```

**Strain-level analysis**

Same here. \

```./02_strain-level_analysis.sh```

**Functional gene content analysis**


