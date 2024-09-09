# GERP

These are helper scripts to assist with Genomic Evolutionary Rate Profiling (GERP).

### brownAlgaeAlignments_splitByScaffold.R

This *R* script takes multiple, separate brown algae genome alignments, splits them by chromosome, and then combines all species' alignments into a single file for each chromosome.

The script is based on but heavily modified from:

https://github.com/BeckySTaylor/Phylogenomic_Analyses/blob/main/GERP_running

Taylor, R.S., Manseau, M., Keobouasone, S., Liu, P., Mastromonaco, G., Solmundson, K., Kelly, A., Larter, N.C., Gamberg, M., Schwantje, H., et al. (2024). High genetic load without purging in caribou, a diverse species at risk. *Curr. Biol.* 34, 1234-1246.e7. https://doi.org/10.1016/j.cub.2024.02.002.

### get_ancestralDerived.R

This *R* script determines which allele is ancestral vs. derived, based on alignments of outgroup reference genomes to the focal species.

### get_derivedAlleleFreqs.R

This *R* script calculates the derived allele frequencies globally, per population, and per individual, based on a VCF of genotypes of all individuals of the focal species.
