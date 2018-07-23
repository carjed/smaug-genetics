## The SMAUG (Singleton Mutation Analysis Utility with Graphics) Pipeline

### Summary

This repository contains the scripts necessary to reproduce the analysis of genome-wide germline mutation patterns described in [Carlson et al., (2017)](http://biorxiv.org/content/early/2017/02/14/108290).

This pipeline can be applied to analyze mutation patterns from any human genomic data mapped to the hg19/GRCh37 reference genome. The only input required is indexed .vcf file(s) containing the standard information (chromosome, position, reference, alternate, etc.) of single-nucleotide variants in a sample. Individual genotypes for each site are not necessary.

### Outline
1. Create a project directory at `/path/to/project/` where we will store this codebase and subsequent output.

1. Enter your project directory and clone this repository with the command: `git clone https://github.com/carjed/smaug-genetics.git`.

1. Edit the `/path/to/project/smaug-genetics/_config.yaml` file to specify the absolute paths to the input directory, project folder, and other parameters.

1. Run `bash download_ref_data.sh` to download the necessary reference data.

1. Run `perl data_mgmt/data_prep/vcf_to_summary.pl copy` to extract positions and alleles of singleton (or common) variants.

1. Run `perl data_mgmt/data_prep/augment_batch.pl` to annotate these summary files with local K-mer sequence context and get data file of motif counts in the reference genome.

1. Run `R/mod_shell.r` to begin analysis for additional processing. This can be run from the command line with `Rscript`, but I recommend running within an R session to troubleshoot any errors without having to reload the data every time (*may have high memory usage \[~30-40GB\]*).

1. The analyses called by `R/mod_shell.r` can be categorized as sequence-only, or sequence+features. The first time you run this script, the sequence-only analyses will run, but if you try to run any of the sequence+features analyses, you will get a notification that you haven't run the proper models. Doing so requires some additional processing, which is handled by scripts in the `data_mgmt` directory.

1. Follow the steps in `data_mgmt/per-site_dp/README.md` to extract depth information from glf files.

1. Follow the steps in `data_mgmt/logit_scripts/README.md` to generate model input data and run the models.

1. Follow the steps in `data_mgmt/process_predicted/README.md` to parse the model output.

1. Assuming the above steps have completed successfully, running `R/mod_shell.r` should execute the final set of analyses.

### Dependencies

#### Hardware
All analyses were run on a server running Ubuntu 14.04. We recommend at least 40GB of RAM and 2-3TB of free disk space for storing output and intermediate file (*disk space recommendation does* NOT *include space used by input files*).

Many of the computationally intensive elements of this pipeline are run using custom batch files processed with the [slurm](http://slurm.schedmd.com/slurm.html) workload management system, which queues and distributes jobs to run simultaneously on multiple server nodes. At peak efficiency, our system utilized 300-400 cores simultaneously, but mileage may vary according to your specific hardware/software configuration, resource availability, etc. For this reason, we do not have a specific recommendation for minimum number of processors...the more, the better!

Note that the slurm batch files generated in this pipeline are included here primarily for reference purposes, but will likely need to be customized to your unique server environment or ported to your workload management system of choice in order to run these steps successfully.

#### Software
- R (v3.2)
  - ggplot2
  - dplyr
  - tidyr
  - broom
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
  - grid
  - gridExtra
  - gtable
- perl (v5.18)
  - POSIX
  - File::Basename
  - File::Path qw(make_path)
  - Benchmark
  - FindBin
  - YAML::XS 'LoadFile'
  - Getopt::Long
  - Math::Round
  - Pod::Usage
  - Compress::Zlib
  - IO::Compress::Gzip
- bash (v4.2)
- [bcftools v1.3.1](http://www.htslib.org/)
- [vcftools](https://vcftools.github.io/index.html)
- [bedtools v2.22.0](http://bedtools.readthedocs.io/en/latest/)
- [samtools-hybrid v0.1.7](https://github.com/statgen/samtools-0.1.7a-hybrid)
  - (this is a modified version of [samtools](http://www.htslib.org/) with better support for .glf files)
- [bigWigToWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig)
- [wigToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig)

We leave it to the user to download these programs and any missing libraries and ensure they run properly. The `R/get_functions.r` script include a helper function named `usePackage()` that automatically installs and loads R packages from the CRAN repository.

### Reference data

This pipeline requires the use of several external data files, almost all of which are publicly available. A script for downloading and formatting these files can be found in this repository at `download_ref_data.sh`. These files include:

- [hg19/GRCh37 reference genome](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz) (all other files must be mapped to this build)
- [chromosome sizes for reference genome](https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)
- RefSeq v69 exon positions
- [Histone modifications (broad peaks) for peripheral blood mononuclear cells from Roadmap Epigenomics Project](http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/)
- [Lamin B1 domains](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz)
- [DNase Hypersensitive Sites](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz)
- [CpG island coordinates](http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt)
  - [Original source (dead link)](http://rafalab.jhsph.edu/CGI/)
- [deCODE sex-averaged recombination rates](http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw)
- [Replication timing profile of lymphoblastoid cell lines](http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip)
- *de novo* mutations from two studies:
  - [Genomes of the Netherlands](https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt)
  - [Inova Translational Medicine Institute Molecular Study of Preterm Birth](http://www.nature.com/ng/journal/v48/n8/extref/ng.3597-S3.xlsx)

### Project directory structure

```
INPUTDIR/
|---vcfs/                       [initial vcf directory]
|---summaries/                  [initial summary files]

PROJDIR/
|---slurm/                      [generated slurm batch files]
|---images/                     [figures output by R scripts]
|---reference_data/             [external files for analysis]
|---smaug-genetics/             [this repository]
    |---data_mgmt/
        |---lib/
        |---data_prep/
        |---logit_scripts/
        |---per-site_dp/
        |---process_predicted/
    |---R/
    |---sandbox/                  [old files]
|---output/
    |---{K}bp_{N}k/               [augmented summary info]
    |---glf_depth/
        |---meandp/
        |---chr*/                   [intermediate depth files]
    |---logmod_data/
        |---coefs/                  [estimated effects of genomic features]
        |---motifs/                 [per-subtype model input]
    |---predicted/                [predicted single-base mutation rates]
        |---tracks/                 [OPTIONAL: predicted rates in .bw format]
    |---slurm/                    [slurm error logs]
```
