This repository contains custom scripts and files for the manuscript:

Bemmels, JB, S Starko, BL Weigel, K Hirabayashi, A Pinch, C Elphinstone, MN Dethier, LH Rieseberg, JE Page, CJ Neufeld, GL Owens. 2025. Population genomics reveals strong impacts of genetic drift without purging and guides conservation of bull and giant kelp. <i>Current Biology</i>, in press.

There was no substantial new program or software developed for this manuscript; rather, the focus of this repository is providing a few data files and several simple scripts to demonstrate calculations, as a supplement to the verbal explanations in the manuscript and manuscript supplementary info.

See the readme files in the subdirectories (contents also pasted below) for further info about the subdirectory contents.

The pages and documents in this repository were developed by Jordan Bemmels (jbemmels@uvic.ca) unless otherwise specified.

# Data

## Brown algae phylogenies

These Newick-formatted phylogenies are for use in Genomic Evolutionary Rate Profiling (GERP) analyses. Two versions are provided that were used in the GERP analyses for bull kelp and giant kelp, respectively.

The phylogenies are composites compiled by hand (with some rescaling, as described in the manuscript) from several sources:

Bringloe, T.T., Starko, S., Wade, R.M., Vieira, C., Kawai, H., De Clerck, O., Cock, J.M., Coelho, S.M., Destombe, C., Valero, M., et al. (2020). Phylogeny and evolution of the brown algae. Crit. Rev. Plant Sci. *39*, 281–321. https://doi.org/10.1080/07352689.2020.1787679

Cánovas, F.G., Mota, C.F., Serrão, E.A., and Pearson, G.A. (2011). Driving south: a multi-gene phylogeny of the brown algal family Fucaceae reveals relationships and recent drivers of a marine radiation. BMC Evol. Biol. *11*, 371. https://doi.org/10.1186/1471-2148-11-371

Choi, S.-W., Graf, L., Choi, J.W., Jo, J., Boo, G.H., Kawai, H., Choi, C.G., Xiao, S., Knoll, A.H., Andersen, R.A., et al. (2024). Ordovician origin and subsequent diversification of the brown algae. Curr. Biol. *34*, 740–754. https://doi.org/10.1016/j.cub.2023.12.069

Kawai, H., Hanyuda, T., Draisma, S.G.A., Wilce, R.T., and Andersen, R.A. (2015). Molecular phylogeny of two unusual brown algae, *Phaeostrophion irregulare* and *Platysiphon glacialis*, proposal of the Stschapoviales ord. nov. and Platysiphonaceae fam. nov., and a re‐examination of divergence times for brown algal orders. J. Phycol. *51*, 918–928. https://doi.org/10.1111/jpy.12332

Silberfeld, T., Leigh, J.W., Verbruggen, H., Cruaud, C., de Reviers, B., and Rousseau, F. (2010). A multi-locus time-calibrated phylogeny of the brown algae (Heterokonta, Ochrophyta, Phaeophyceae): investigating the evolutionary nature of the “brown algal crown radiation.” Mol. Phylogenet. Evol. *56*, 659–674. https://doi.org/10.1016/j.ympev.2010.04.020

Starko, S., Soto Gomez, M., Darby, H., Demes, K.W., Kawai, H., Yotsukura, N., Lindstrom, S.C., Keeling, P.J., Graham, S.W., and Martone, P.T. (2019). A comprehensive kelp phylogeny sheds light on the evolution of an ecosystem. Mol. Phylogenet. Evol. *136*, 138–150. https://doi.org/10.1016/j.ympev.2019.04.012

## Coastline rasters

These rasters were used for calculating the geographic distance 'as the kelp floats' (i.e., the least-cost distance travelled by ocean) between kelp populations. Each pixel in the the southern BC and Washington raster represents one millidegree, whereas in the northern BC raster each pixel represents two millidegrees.

Note that some small areas of the rasters have been modified by hand to convert non-ocean pixels to ocean pixels. The purpose of this was to ensure that all kelp populations occurred in ocean pixels, and to ensure that narrow ocean passages were fully open to dispersal by kelp.

Rasters were generated from a polygon of the BC and partial Washington coastline:

GeoBranch BC. 2002. NTS BC Coastline Polygons 1:250,000 - Digital Baseline Mapping (NTS). *British Columbia Data Catalogue.* Downloaded from [https://catalogue.data.gov.bc.ca/dataset/nts-bc-coastline-polygons-1-250-000-digital-baseline-mapping-nts](https://catalogue.data.gov.bc.ca/dataset/nts-bc-coastline-polygons-1-250-000-digital-baseline-mapping-nts) \[accessed 2022/06/21\].

A derivative product of the original polygons is being reposted here in accordance with an [Open Government License - Canada v.2.0](https://open.canada.ca/en/open-government-licence-canada).

# Scripts

## GERP

These are helper scripts to assist with Genomic Evolutionary Rate Profiling (GERP).

### brownAlgaeAlignments_splitByScaffold.R

This *R* script takes multiple, separate brown algae genome alignments, splits them by chromosome, and then combines all species' alignments into a single file for each chromosome.

The script is based on but heavily modified from:

https://github.com/BeckySTaylor/Phylogenomic_Analyses/blob/main/GERP_running

Taylor, R.S., Manseau, M., Keobouasone, S., Liu, P., Mastromonaco, G., Solmundson, K., Kelly, A., Larter, N.C., Gamberg, M., Schwantje, H., et al. (2024). High genetic load without purging in caribou, a diverse species at risk. *Curr. Biol.* 34, 1234-1246.e7. https://doi.org/10.1016/j.cub.2024.02.002.

### get_ancestralDerived.R

This *R* script determines which allele is ancestral vs. derived, based on alignments of outgroup reference genomes to the focal species.

## SnpEff

These are helper scripts to assist with SnpEff analysis.

### vcf_to_alleleFreq.sh

This bash script will process a VCF file that has been annotated by SnpEff and print a custom allele frequency text file that has the frequencies of reference and alternative alleles for each site, plus the SnpEff annotation.

### add_ancestralDerived_to_SnpEff.R

This *R* script will add the ancestral and derived definitions (previously defined during GERP analysis) to the custom allele frequency text file.

### getDerivedAlleleFreqs_snpEff

This *R* script calculates the derived allele frequencies globally, per population, and per individual, based on a VCF of genotypes of all individuals of the focal species.

## Purging and drift

These scripts calculate allele frequencies by population in different site categories, which are needed to statistically test the predictions of purging and genetic drift.

### purging_drift_GERP.R

This *R* script calculates the allele frequencies for evolutionarily labile and evolutionarily conserved sites, from the GERP analysis.

### purging_drift_SnpEff.R

This *R* script calculates the allele frequencies for modifier and low-, moderate-, and high-impact sites, from the SnpEff analysis.

## Realized load of empirical individuals

These scripts calculate realized genetic load of empirically observed individuals.

### realizedLoad_empirical_GERP.R

This *R* script calculates the empirical realized load for evolutionarily labile and evolutionarily conserved sites, from the GERP analysis.

### realizedLoad_empirical_SnpEff.R

This *R* script calculates the empirical realized load for modifier and low-, moderate-, and high-impact sites, from the SnpEff analysis.

## Realized load of simulated crosses

These scripts calculate realized genetic load of individuals generated from simulated crosses between and within populations.

### getGTs_ancDer_GERP.R

This *R* script prints and pre-processes the genotypes for all individuals in each population, converting genotypes to [0, 1, 2] format indicating the number of copies of the derived allele, in preparation for the GERP analysis.

### realizedLoad_simulatedCrosses_GERP.R

This *R* script calculates the realized load from simulated crosses for evolutionarily labile and evolutionarily conserved sites, from the GERP analysis. 

### getGTs_ancDer_SnpEff.R

This *R* script prints and pre-processes the genotypes for all individuals in each population, converting genotypes to [0, 1, 2] format indicating the number of copies of the derived allele, in preparation for the SnpEff analysis.

### realizedLoad_simulatedCrosses_SnpEff.R

This *R* script calculates the realized load from simulated crosses for modifier and low-, moderate-, and high-impact sites, from the SnpEff analysis. 

## Runs of Heterozygosity (ROHets)

These scripts are to identify Runs of Heterozygosity (ROHets), not to be confused with Runs of Homozygosity (ROHs). ROHets are short regions of heterozygosity surrounded by ROHs. Identifying short ROHets and masking windows of the genome that are exceptionally frequently found in an ROHet was done prior to estimating effective populations size from ROHs, under the assumption that windows that are frequently found in ROHets may be confounded regions of the genome or represent errors. This would be problematic as unmasked ROHet regions could interfere with the correct inference of the size distribution of ROHs.

### identify_short_ROHets.R

This *R* script identifies regions of the genome that are in ROHets for each population.

### mask_ROHets.R

This *R* script creates a mask to exclude windows that are found in ROHets more frequently (across all individuals) than the 99.99th percentile expectation. The mask can then be used to re-run *bcftools roh*.


### get_derivedAlleleFreqs.R

This *R* script calculates the derived allele frequencies globally, per population, and per individual, based on a VCF of genotypes of all individuals of the focal species.

