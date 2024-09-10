# Realized load of simulated crosses

These scripts calculate realized genetic load of individuals generated from simulated crosses between and within populations.

### getGTs_ancDer_GERP.R

This *R* script prints and pre-processes the genotypes for all individuals in each population, converting genotypes to [0, 1, 2] format indicating the number of copies of the derived allele, in preparation for the GERP analysis.

### realizedLoad_simulatedCrosses_GERP.R

This *R* script calculates the realized load from simulated crosses for evolutionarily labile and evolutionarily conserved sites, from the GERP analysis. 

### getGTs_ancDer_SnpEff.R

This *R* script prints and pre-processes the genotypes for all individuals in each population, converting genotypes to [0, 1, 2] format indicating the number of copies of the derived allele, in preparation for the SnpEff analysis.

### realizedLoad_simulatedCrosses_SnpEff.R

This *R* script calculates the realized load from simulated crosses for modifier and low-, moderate-, and high-impact sites, from the SnpEff analysis. 
