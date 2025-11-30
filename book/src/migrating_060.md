# Migrating to Modkit v0.6.0 `pileup`.

## Changes to pileup algorithm

1. `--ignore` option has been removed. This option was often a source of confusion and has been removed in favor of the following alternatives:
  - Use `--modified-bases` to emit records for only the base modifications you're interested in.
  For example, a direct RNA run may have many base modifications on each of the primary sequence bases.
  If you're only interested in m6A use `--modified-bases m6A` and only get records with this base modification.
  - Use `--combine-mods` when you want to know if a base is "modified vs unmodified".
  This can be necessary when comparing to results where all base modifications may be marked as a single base modification.
1. Use `--phased` instead of `--partition-tags`. 
  The `--phased` flag is a more performant version of `--partition-tags HP`.
1. The `--bedgraph` flag has been removed.
  The recommended way to create a "trace"-like output file is to make a bedMethyl and either pipe it directly into `modkit bedmethyl tobigwig` or run the `tobigwig` command following bedMethyl creation.
1. Duplicate read names are no longer logged.
1. `--max-depth` removed. Modkit v0.6.0 has a new algorithm that doesn't use `max-depth`.
If you have a modBAM with **very high** depth and you don't want to tabulate counts for this depth it is currently recommended to subsample the modBAM before using Modkit.
1. `--edge-filter`-ed base modification counts are no longer tabulated in N<sub>nocall</sub>, they will be completely ignored from tabulation.
1. Simpler threshold mechanism when using `--modified-bases`.
  Prior to v0.6.0 the threshold evaluation was to take the most likely passing base modification call.
  This lead to confusion about how to "tune" this parameter to tweak results.
  When using `--modified-bases` (strongly recommended) the algorithm is now simply to use the most likely call (the one with the highest probability) and determine if it is pass or fail by whether or not it is less than the threshold for that base or modification.
  In almost all cases this simpler algorithm produces the same results and is what the user expected in the first place.
1. `--combine-strands` when using multiple palindromic motifs is no longer implemented.
  If this is a crucial step in your workflow, please open a GitHub issue.
  The new workflow is to run with a single motif and `--combine-strands` for each motif.
  Modkit v0.6.0 is substantially quicker and lighter on resource requirements such that running multiple times will likely be quicker than the older versions.
1. `--invert-edge-filter` removed.
  This option is still available in `modkit modbam adjust-mods`.

