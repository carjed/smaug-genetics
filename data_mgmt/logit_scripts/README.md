# Modeling effect of sequence context and genomic features

Our model to predict the effects of sequence context and genomic features requires a series of preprocessing steps to generate a set of 24,576 files (1 per each 7-mer subtype), to be passed to `R/log_mod.r` in our analyses.

## Preparing input data
----------------------------------------------

*Before running the scripts in this directory, we must have already generated indices for per-site average depth, which should be located at* `output/glf_depth/meandp/chr*.dp`. *See `data_mgmt/per-site_dp/README.md` for details on how to generate these files.*

*Also, when running `R/mod_shell.r` the first time, make sure the `build_logit` flag in `_config.yaml` is set to* TRUE *to index the singleton sites when building this model. These files will be located at `output/logmod_data/chr*_{type}_sites.txt`.*

----------------------------------------------

## Running the scripts

Once the necessary depth and singleton indices have been generated, we will run scripts to execute three batch files:

1. `build_data_batch.pl`
2. `split_data_batch.pl`
3. `runmod_batch.pl`

### Building data

The `build_data_batch.pl` script will generate and execute a slurm batch file, which writes a table for each chromosome/type with 5 columns:
1. CHR
2. POS
3. SEQUENCE \[7-mer motif centered at site, with symmetric sequence\]
4. MUT \[indicates if site is singleton (1) or non-mutated (0)\]
5. DP \[average depth of coverage at site\]

These compressed files (~55GB total) are located at `output/logmod_data/chr*_${type}_full.txt.gz`.

### Splitting by subtype

Because our model runs independently for each subtype, we need to split the input data by type/motif by running `split_data_batch.pl`. This creates and executes a slurm batch file that writes subtype-specific input data to `/output/logmod_data/motifs/{type}_{motif}.txt`.

Note that unlike the `build_data_batch.pl` script, `split_data_batch.pl` executes a single batch file with 6 jobs in the array (one per type, each sequentially looping through 22 chromosomes) rather than 6 batch files with 22 jobs in the array that run simultaneously. This is because when each `chr*_${type}_full.txt.gz` file is split, the data for a given subtype is appended to the same `{type}_{motif}.txt` file. Consequently, this script can take a long time to run.

## Running the models
Once the subtype-specific input files have been generated, we are ready to run the models. The `runmod_batch.pl` script creates and executes 6 slurm batch files (one per type), calling `R/log_mod.r` independently for each of the 24,576 7-mer subtypes.

For each subtype, we obtain two output files:
1. the beta coefficients and standard errors for the features considered in the model, located at `output/logmod_data/coefs/{type}/{type}_{motif}_coefs.txt`
2. the predicted mutation rates for all sites of that subtype. We will eventually need these rates to be ordered by chromosome/position, so the predicted rates are split by chromosome into `output/predicted/{type}/chr*/{type}_{motif}.txt`. (See `data_mgmt/process_predicted/README.md` for details on how we recombine these into per-chromosome/per-type files).

Because `runmod_batch.pl` submits several thousand jobs to slurm, there is a strong possibility that many will to fail due to server load, queue limits, etc. Running this script with the `--parentjob [job number]` option will scan the slurm monitor for failed jobs and generate/execute secondary batch files to re-run them. Note that you may need to manually adjust the parameters in the batch file if they were the cause of failure (e.g., exclude problematic nodes, set memory limits).
