##### Script is to print the genotypes for all individuals in each population. Genotypes are recoded in [0, 1, 2] format relative to ancestral-derived definitions, indicating whether there are 0, 1, or 2 copies of the derived allele. This script processes the SnpEff analysis. #####

# Author use only: GitHub version derived from version 5 of KL14 SnpEff script 01b_getGTs_KL14_KRI_ancestralGERP_unprocessedPopsOnly.R, 2024/09/10;

# This example is to determine the derived allele frequencies for bull kelp (Nereocystis). If the focal species in instead giant kelp (Macrocystis), the user needs to change "Macrpy" to "Nerelu" throughout this script, so that the outgroup will be set to Nereocystis luetkeana.;

#####

### Prerequisites:

# 1) List of all samples, one individual per line;

	# sampleList.txt;

# 2) Printed genotypes for all samples in the VCF, this can be obtained as follows:

	# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -R sitesToPrint.txt file.bcf > genotypes.txt

# 3) Printed genotypes are split into smaller files of 1,000 lines each, for efficient processing:

	# split -l 1000 genotypes.txt genotypes_splitFile/geno_

# 4) Allele-frequency file for ancestral and derived alleles (previously added with the getDerivedAlleleFreqs_snpEff.R example script) are split into smaller files of 1,000 lines each, with exactly the same sites in the same order in each of the smaller files (tail +2 avoids printing the header, to ensure that sites are identical):

	# tail +2 alleleFreq_ancestralDerived.txt | split -l 1000 - ancestralDerived_splitFile/anc_
		
##############################
##############################

system("mkdir tmp", inter = F);
system("mkdir genotypes_snpEff_byPop", intern = F);

# define all samples and their populations;

samples <- read.table("sampleList.txt", header = F, sep = "\t", stringsAsFactors = F);
colnames(samples) <- "ID";

# user may need to edit these lines if the individuals and populations are encoded in the sample names in a different way;
samples$individual <- gsub("KL..\\.", "", samples$ID);
samples$pop <- substr(samples$individual, 1, 8);

# summarize the information about number of samples per population;

samples_byPop <- data.frame(pop = unique(samples$pop), n_total = NA);

for (i in 1:nrow(samples_byPop)) {
	samples_byPop$n_total[i] = sum(samples$pop == samples_byPop$pop[i]);
}

# define the populations to process - this would be any populations that have at least three individuals present;

allPops_min3ind <- unique(samples_byPop$pop[samples_byPop$n_total >= 3]);

##############################
##############################

library(foreach);
library(doParallel);

##### Process one split-file at a time #####

allGenoFiles <- list.files("genotypes_splitFile/", pattern = "geno_");
allSplits <- gsub("geno_", "", allGenoFiles);

#

registerDoParallel(cores=length(allPops_min3ind));

foreach (k=1:length(allPops_min3ind)) %dopar% {
  
  cur_pop <- allPops_min3ind[k];
  
  # we want to print the genotypes for bothe selfed and non-selfed individuals;
  cur_samples <- samples$individual[samples$pop == cur_pop];
  
  print(paste0("Working on population ", cur_pop));
  
  for (i in 1:length(allSplits)) {
    
    print(paste0("Working on split file ", i, " of ", length(allSplits)));
    
    cur_anc <- read.table(paste0("ancestralDerived_splitFile/anc_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
    colnames(cur_anc) <- c("chr", "pos", "ref", "alt", "ancestral", "derived", "freqAlt", "freqDer", "annotation", "putative_impact", "gene_name");
    
    cur_geno <- read.table(paste0("genotypes_splitFile/geno_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
    colnames(cur_geno) <- c("chr", "pos", "ref", "alt", samples$individual);
    
    # subset to only the genotypes for individuals in the current population;
    cur_geno <- cur_geno[ , c("chr", "pos", "ref", "alt", cur_samples)];
    
    cur_df <- merge(cur_anc[ , c("chr", "pos", "ref", "alt", "ancestral", "derived", "annotation", "putative_impact")], cur_geno, by = c("chr", "pos", "ref", "alt"), sort = F);
    
    # change from 0/0, 0/1, 1/1 format to 0, 1, 2 format, where 0 = homozygous for ANCESTRAL allele, 1 = heterozygous, 2 = homozygous for DERIVED allele;
    cur_df$ref_is_anc <- NA;
    cur_df$ref_is_anc[cur_df$ref != cur_df$anc] <- 0;
    cur_df$ref_is_anc[cur_df$ref == cur_df$anc] <- 1;
    #
    cur_df[cur_df$ref_is_anc == T & cur_df == "0/0"] <- 0;
    cur_df[cur_df$ref_is_anc == T & cur_df == "0/1"] <- 1;
    cur_df[cur_df$ref_is_anc == T & cur_df == "1/1"] <- 2;
    cur_df[cur_df$ref_is_anc == T & cur_df == "./."] <- NA;
    #
    cur_df[cur_df$ref_is_anc == F & cur_df == "0/0"] <- 2;
    cur_df[cur_df$ref_is_anc == F & cur_df == "0/1"] <- 1;
    cur_df[cur_df$ref_is_anc == F & cur_df == "1/1"] <- 0;
    cur_df[cur_df$ref_is_anc == F & cur_df == "./."] <- NA;
    #
    # we have some sites that although the outgroup species has a reference allele, the ancestral allele could not be defined because the outgroup's allele is different from either the ref or the alt allele;
    # anything left by now should be changed to NA;
    cur_df[cur_df == "0/0"] <- NA;
    cur_df[cur_df == "0/1"] <- NA;
    cur_df[cur_df == "1/1"] <- NA;
    cur_df[cur_df == "./."] <- NA;
    #
    write.table(cur_df[ , c("chr", "pos", "annotation", "putative_impact", cur_samples)], paste0("tmp/", cur_pop, "_", allSplits[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }
  
  # combine into a single file;
  # do not print the first line, as it is a header in all files - print the header first as its own separate line;
  # then, in the second command, append using ">>";
  
  system(paste0("cat tmp/", cur_pop, "_aa.txt | head -1 > genotypes_snpEff_byPop/genotypes_derAnc_snpEff_", cur_pop, ".txt"));
  system(paste0("tail -n +2 -q tmp/", cur_pop, "* >> genotypes_snpEff_byPop/genotypes_derAnc_snpEff_", cur_pop, ".txt"));
  
  # delete the temporary files;
  system(paste0("rm tmp/", cur_pop, "*"));
  
}
