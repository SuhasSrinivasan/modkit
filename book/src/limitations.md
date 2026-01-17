# Current limitations

Known limitations and forecasts for when they will be removed.

1. During `modkit pileup`, it is assumed that each read should only have one primary alignment.
  If a read has more than one primary alignment it will be double counted, it is therefore recommended to remove duplicates before running Modkit.
1. The MAP-based p-value metric ([details](./dmr_scoring_details.md#map-based-p-value)) performs a test that there is a difference in modification (of any kind) between two conditions.
  If a position has multiple base modification calls (such as 5hmC and 5mC) the calls are summed together into a single "modified" count. 
  If a position differs only in the modification _type_ (such as one condition has more 5hmC and the other has more 5mC) this effect will not be captured in the MAP-based p-value.
  The [likelihood ratio](./dmr_scoring_details.md#likelihood-ratio-scoring-details) test _does_ capture changes in modification type.
1. The MAP-based p-value is not available when performing DMR on regions.
  This is because there is potentially large variability in the number of modified bases and coverage over regions.
  This variability translates into varying degrees of statistical power and makes comparisons difficult to interpret.
1. As of v0.6.0 `pileup` will saturate depth at 65,535 (maximum for a unsigned 16-bit integer). If your modBAM has a depth greater than this value, it is recommended to use the `--high-depth` flag so that 65,535 reads will be used at each genomic position. To achieve greater depth, divide the modBAM into subsets with depth <= 65,535, run `pileup` on each, then `modbam merge`. This limitation will be removed in the future when it can be proven not to have any performance or resource regression.
