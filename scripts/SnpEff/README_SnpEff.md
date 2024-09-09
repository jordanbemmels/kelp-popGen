# SnpEff

These are helper scripts to assist with SnpEff analysis.

### vcf_to_alleleFreq.sh

This bash script will process a VCF file that has been annotated by SnpEff and print a custom allele frequency text file that has the frequencies of reference and alternative alleles for each site, plus the SnpEff annotation.

### add_ancestralDerived_to_SnpEff.R

This *R* script will add the ancestral and derived definitions (previously defined during GERP analysis) to the custom allele frequency text file.

### getDerivedAlleleFreqs_snpEff

This *R* script calculates the derived allele frequencies globally, per population, and per individual, based on a VCF of genotypes of all individuals of the focal species.
