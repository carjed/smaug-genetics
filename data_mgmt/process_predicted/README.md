# Processing per-site predicted rates

## Merging and sorting predicted rates
Because the predicted rates are written as per-chromsome, per-subtype files, we must concatenate all subtype files and sort by position for each chromsome. This is performed by running `sort_pred_batch.pl`, which creates and executes 6 slurm batch files (one per type), generating per-chromosome/per-type files into `output/predicted/chr*.{type}.txt`.

## Extracting null background for validation models
These predicted rate files are the basis for generating the non-mutated background we use in our validation models on *de novo* mutations. The `sample_mu.pl` script is responsible for subsampling ~30 million non-mutated sites (~1% of the genome) and merges these sites with the *de novo* sites into a file named `output/validation_sites.txt` to be passed to the `validation.r` script for further analysis.

## Creating UCSC Genome Browser-compatible tracks (optional)
With an additional processing step, we can convert the predicted rate tables into BigWig (.bw) format for compatibility with the UCSC Genome Browser. This is performed using the `toWig.pl` script (if desired, a slurm batch file to do this across all chromosomes is available at  `data_mgmt/toWig/toWig_batch.txt`).

## Binning predicted rates (optional)
We may be interested in the total expected number of singletons in a genomic window. This is done by running `bin_pred_batch.pl`, which creates and executes 6 batch files to obtain per-bin expected counts for each type, located in `output/predicted/binned/chrs/chr*_{type}_binned.txt`.
