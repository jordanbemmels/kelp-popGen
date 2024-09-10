# Runs of Heterozygosity (ROHets)

These scripts are to identify Runs of Heterozygosity (ROHets), not to be confused with Runs of Homozygosity (ROHs). ROHets are short regions of heterozygosity surrounded by ROHs. Identifying short ROHets and masking windows of the genome that are exceptionally frequently found in an ROHet was done prior to estimating effective populations size from ROHs, under the assumption that windows that are frequently found in ROHets may be confounded regions of the genome or represent errors. This would be problematic as unmasked ROHet regions could interfere with the correct inference of the size distribution of ROHs.

### identify_short_ROHets.R

This *R* script identifies regions of the genome that are in ROHets for each population.

### mask_ROHets.R

This *R* script creates a mask to exclude windows that are found in ROHets more frequently (across all individuals) than the 99.99th percentile expectation. The mask can then be used to re-run *bcftools roh*.
