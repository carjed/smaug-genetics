## The SMAUG (Singleton Mutation Analysis Utility with Graphics) Pipeline

### Summary

This repository contains the scripts necessary to reproduce the analysis of genome-wide germline mutation patterns described in (citation here). Scripts are written in R, perl, and UNIX shell.

### Dependencies

#### Hardware
All analyses were run on a server running Ubuntu 14.04. We recommend at least 5TB of free disk space for storing input/output/and intermediate files, and 40GB of RAM. Many of the computationally intensive elements of this pipeline are run using custom batch files processed with the [slurm](http://slurm.schedmd.com/slurm.html) workload management system, which queues and distributes jobs to run simultaneously on different nodes of our server cluster. At peak efficiency, our system utilized 300-400 cores simultaneously, but mileage may vary according to your specific hardware/software configuration, resource availability, etc. For this reason, we do not have a specific recommendation for minimum number of processors...the more, the better!

The original slurm batch files used in this pipeline (and any scripts used to generate them) are included here for reference purposes, but must be customized to your unique server environment or implemented in an alternative system in order to run these steps successfully.

#### Software
- R (version 3.2)
  - ggplot2
  - dplyr
  - tidyr
  - reshape2
  - RColorBrewer
  - MASS
  - speedglm
  - boot
  - devtools
  - psych
  - lmtest
  - fmsb
  - cowplot
  - hexbin
  - stringr
- perl (version 5.14)
- bash (version 4.2)
- [bcftools](http://www.htslib.org/)
- [vcftools](https://vcftools.github.io/index.html)
- [bedtools](http://bedtools.readthedocs.io/en/latest/)
- [samtools-hybrid](https://github.com/statgen/samtools-0.1.7a-hybrid)
  - (this is a modified version of [samtools](http://www.htslib.org/) with better support for .glf files)
- [bigWigToWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig)
- [wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig)

We leave it to the user to download these programs and ensure they run properly.


### Reference data

This pipeline requires the use of several external data files, all publicly available. A script for downloading and formatting these files can be found in this repository at `download_ref_data.sh`. These files include:

- [hg19/GRCh37 reference genome](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz) (all other files must be mapped to this build)
- [chromosome sizes for reference genome](https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)
- RefSeq exon positions
- [Histone modifications (broad peaks) for peripheral blood mononuclear cells from Roadmap Epigenomics Project](http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/)
- [Lamin B1 domains](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz)
- [DNase Hypersensitive Sites](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz)
- [CpG island coordinates](http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt)
  - [Original source (dead link)](http://rafalab.jhsph.edu/CGI/)
- [deCODE sex-averaged recombination rates](http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw)
- [Replication timing profile of lymphoblastoid cell lines](http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip)
- *de novo* mutations from two studies:
  - Genomes of the Netherlands
  - Inova Newborn Screening Cohort

### Data preparation

This pipeline can be applied to estimate mutational spectra from any human genomic data mapped to the hg19/GRCh37 reference genome. The only input required is an indexed .vcf file containing the standard information (chromosome, position, reference, alternate, etc.) of single-nucleotide variants in a sample. Individual genotypes for each site are not necessary.

### Outline

1. Enter your project directory and clone this repository

1. Edit the `options.txt` file to specify

1. Download the necessary reference data:
  `bash PROJDIR/data_mgmt/process_covs/download_ref_data.sh`

1. Start with an input VCF, and get per-chromosome summary of singleton (or common) variants using bcftools (data_pipeline/vcf_to_summary.pl)

1. Annotate these summary files with local K-mer sequence context, count  K-mers across the genome (data_pipeline/augment_summary.pl). *This script is implemented in a slurm batch file for more efficient utilization of server resources.

1. These annotated summaries are passed to an R script (R/mod_shell.r) for additional processing.
-Get summary of K-mer relative mutation rates
-Estimate effects of genomic features and predict position-specific mutation rates using logistic regression

### Project directory structure

```
INPUTDIR\
|---vcfs\
|---summaries\

PROJDIR\
|---smaug-genetics\ [scripts from this repository]
    |---data_mgmt\
        |---logit_scripts\
        |---per-site_dp\
        |---process_covs\
        |---process_dnms\
    |---data_pipeline\
    |---R\
|---images\ [contains all figures output by R scripts]
|---reference_data\ [contains all reference data downloaded by download_ref_data.sh script]
    |---histone_marks\
        |---broad\
|---output\
    |---{K}bp_{N}k\ [contains augmented summary info for specified motif/bin lengths]
    |---glf_depth\
    |---logmod_data\
        |---chr*\ [temporary files]
        |---coefs\
        |---motifs\
    |---predicted\
        |---chr*.{type}.txt [predicted single-base mutation rates]
        |---tracks\
```
