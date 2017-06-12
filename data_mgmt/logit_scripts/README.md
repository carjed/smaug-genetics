# Modeling effect of sequence context and genomic features

Our model to predict the effects of sequence context and genomic features requires a series of preprocessing steps to generate a set of 24,576 files (1 per each 7-mer subtype), to be passed to `R/log_mod.r` in our analyses.

## Preparing input data
----------------------------------------------

*Before running the scripts in this directory, we must have already generated indices for per-site average depth, which should be located at* `output/glf_depth/meandp/chr*.txt`. *See `data_mgmt/per-site_dp/README.md` for details on how to generate these files.*

*Also, when running `R/mod_shell.r` the first time, make sure the `build_logit` flag in `_config.yaml` is set to* TRUE *to index the singleton sites when building this model. These files will be located at `output/logmod_data/chr*_sites.txt`.*

----------------------------------------------

## Running the scripts

### Generating input data
Once the necessary depth and singleton indices have been generated, run the `build_data_worker.pl` script to prepare the input data for the models.

This script writes a file for each 5Mb chunk of sequence containing 10 columns:
1. CHR
2. POS
3. SEQUENCE \[7-mer motif centered at site, with symmetric sequence\]
4. AT_CG \[indicator if site is an AT>CG singleton (1) or non-mutated (0)\]
5. AT_GC
6. AT_TA
7. GC_AT
8. GC_CG
9. GC_TA
5. DP \[average depth of coverage at site\]

These files are written to: `output/logmod_data/chr*/chr*.{start}-{end}.txt`. Because our model eventually needs to run independently for each subtype, after each file has been generated, it is divided by the SEQUENCE column and appended into motif-specific files at `output/logmod_data/motifs/{motif}.txt`.

Runtime for this script is ~12-13 hours, though this may be bottlenecked significantly by your available RAM (for the sort operation), disk read/write speeds, and single-core CPU speed.

### Running the models
Once the subtype-specific input files have been generated, we are ready to run the models. The `runmod_batch.pl` script creates and executes 6 slurm batch files (one per type), calling `R/log_mod.r` independently for each set of 4,096 7-mer subtypes corresponding to a given type.

For each subtype, we obtain two output files:
1. the beta coefficients and standard errors for the features considered in the model, located at `output/logmod_data/coefs/{type}/{type}_{motif}_coefs.txt`
2. the predicted mutation rates for all sites of that subtype. We will eventually need these rates to be ordered by chromosome/position, so the predicted rates are split by chromosome into `output/predicted/{type}/chr*/{type}_{motif}.txt`. (See `data_mgmt/process_predicted/README.md` for details on how we recombine these into per-chromosome/per-type files).

Because `runmod_batch.pl` submits several thousand jobs to slurm, there is a strong possibility that many will to fail due to server load, queue limits, etc. Running this script with the `--parentjob [job number]` option will scan the slurm monitor for failed jobs and generate/execute secondary batch files to re-run them. Note that you may need to manually adjust the parameters in the batch file if they were the cause of failure (e.g., exclude problematic nodes, set memory limits).
