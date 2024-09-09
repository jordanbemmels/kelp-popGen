##### This script will calculate the derived allele frequencies at each site (a) globally, (b) for each population, and (c) for each individual. The derived allele frequencies are calculated from called genotypes in a VCF file. The script is run after the get_ancestalDerived.R script has been used to identify which allele is ancetral and which is derived. There are several prerequisites (that also show the input file organization) described below. #####

# Author use only: GitHub version derived from KL14 script getDerivedAlleleFreqs_KL14_ancestralGERP.R, 2024/09/09;

# This example is to determine the derived allele frequencies for bull kelp (Nereocystis). If the focal species in instead giant kelp (Macrocystis), the user needs to change "Macrpy" to "Nerelu" throughout this script, so that the outgroup will be set to Nereocystis luetkeana.;

### Prerequisites:

# 1) List of all samples (individuals) in the VCF file, this can be obtained as follows:

	# bcftools query -l file.bcf > sampleList.txt

# 2) Printed genotypes for all samples in the VCF, this can be obtained as follows:

	# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -R sitesToPrint.txt file.bcf > genotypes.txt

# 3) Printed genotypes are split into smaller files of 1,000 lines each, for efficient processing:

	# split -l 1000 genotypes.txt genotypes_splitFile/geno_

# 4) Ancestral and derived alleles (determined previously with get_ancestralDerived.R example script) are split into smaller files of 1,000 lines each, with exactly the same sites in the same order in each of the smaller files (tail +2 avoids printing the header, to ensure that sites are identical):

	# tail +2 ancestralDerived.txt | split -l 1000 - ancestralDerived_splitFile/anc_

##############################
##############################

#####

system(paste0("mkdir tmp"), intern = F);

# define all samples and their populations;

samples <- read.table("sampleList.txt", header = F, sep = "\t", stringsAsFactors = F);
colnames(samples) <- "ID";

# identify the individual and the population based on the sample name;
# NOTE: this should be edited by the user if the user uses a different naming scheme for samples;
samples$individual <- gsub("KL..\\.", "", samples$ID);
samples$pop <- substr(samples$individual, 1, 8);

allPops <- unique(samples$pop);
allInd <- unique(samples$individual);

##############################
##############################

library(foreach);
library(doParallel);

##### Analyze global allele frequencies across all populations #####

##### Process one split-file at a time #####

allGenoFiles <- list.files("genotypes_splitFile/", pattern = "geno_");
allSplits <- gsub("geno_", "", allGenoFiles);

# set up 48 to be run at a time;

n_batches <- 48;
length_batches <- ceiling(length(allSplits) / n_batches);
for (i in 1:n_batches) {
  if (i == 1) {
    batches = list((length_batches*(i-1) + 1):min(length(allSplits), length_batches * i));
  } else {
    batches <- append(batches, list((length_batches*(i-1) + 1):min(length(allSplits), length_batches * i)));
  }
}

#

registerDoParallel(cores=48);

for (k in 1:length(batches)) {
  
  foreach (i=batches[[k]]) %dopar% {
    
    print(paste0("Working on split file ", i, " of ", length(allSplits)));
    
    cur_anc <- read.table(paste0("ancestralDerived_splitFile/anc_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
    colnames(cur_anc) <- c("chr", "pos", "ref", "alt", "N", "S", "Macrpy", "Saccja", "Lamidi", "n_outgroup_alleles", "unique_outgroup_allele", "ancestral", "derived");
    
    cur_geno <- read.table(paste0("genotypes_splitFile/geno_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
    colnames(cur_geno) <- c("chr", "pos", "ref", "alt", samples$individual);
    
    cur_df <- merge(cur_anc[ , c("chr", "pos", "ref", "alt", "N", "S", "ancestral", "derived")], cur_geno, by = c("chr", "pos", "ref", "alt"), sort = F);
    
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
    cur_df$n_homAnc <- apply(cur_df[ , samples$individual], 1, function(x) {sum(x == 0, na.rm = T)});
    cur_df$n_het <- apply(cur_df[ , samples$individual], 1, function(x) {sum(x == 1, na.rm = T)});
    cur_df$n_homDer <- apply(cur_df[ , samples$individual], 1, function(x) {sum(x == 2, na.rm = T)});
    cur_df$n_total <- cur_df$n_homAnc + cur_df$n_het + cur_df$n_homDer;
    #
    # to save file space, we can calculate this later when we read in the file, instead of printing it now (also has long decimals, waste of file space);
    #cur_df$n_ind <- cur_df$n_total;
    #cur_df$freq_ancestral <- round((2*cur_df$n_homAnc + cur_df$n_het) / (2*cur_df$n_total), 4);
    #cur_df$freq_derived <- round((cur_df$n_het + 2*cur_df$n_homDer) / (2*cur_df$n_total), 4);
    #
    write.table(cur_df[ , c("chr", "pos", "N", "S", "n_homAnc", "n_het", "n_homDer")], paste0("tmp/freqs_", allSplits[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
  
  }

}

# combine into a single file;
# do not print the first line, as it is a header in all files - print the header first as its own separate line;
# then, in the second command, append using ">>";

system(paste0("cat tmp/freqs_aa.txt | head -1 > derivedAlleleFreqs.txt"));
system(paste0("tail -n +2 -q tmp/freqs* >> derivedAlleleFreqs.txt"));

# delete the temporary files;
system(paste0("rm tmp/freqs_*"));

##### DONE, IT'S WORKING #####

#####
#####
#####

##### REPEAT FOR EACH POPULATION #####

##### Process one split-file at a time #####

allGenoFiles <- list.files("genotypes_splitFile/", pattern = "geno_");
allSplits <- gsub("geno_", "", allGenoFiles);

#

registerDoParallel(cores=length(allPops));

foreach (k=1:length(allPops)) %dopar% {
  
  cur_pop <- allPops[k];
  cur_samples <- samples$individual[samples$pop == cur_pop];
  
  print(paste0("Working on population ", cur_pop));
  
  for (i in 1:length(allSplits)) {
    
    print(paste0("Working on split file ", i, " of ", length(allSplits)));
    
    cur_anc <- read.table(paste0("ancestralDerived_splitFile/anc_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
    colnames(cur_anc) <- c("chr", "pos", "ref", "alt", "N", "S", "Macrpy", "Saccja", "Lamidi", "n_outgroup_alleles", "unique_outgroup_allele", "ancestral", "derived");
    
    cur_geno <- read.table(paste0("genotypes_splitFile/geno_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
    colnames(cur_geno) <- c("chr", "pos", "ref", "alt", samples$individual);

    # subset to only the genotypes for individuals in the current population;
    cur_geno <- cur_geno[ , c("chr", "pos", "ref", "alt", cur_samples)];
    
    cur_df <- merge(cur_anc[ , c("chr", "pos", "ref", "alt", "N", "S", "ancestral", "derived")], cur_geno, by = c("chr", "pos", "ref", "alt"), sort = F);
    
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
    cur_df$n_homAnc <- apply(cur_df[ , cur_samples], 1, function(x) {sum(x == 0, na.rm = T)});
    cur_df$n_het <- apply(cur_df[ , cur_samples], 1, function(x) {sum(x == 1, na.rm = T)});
    cur_df$n_homDer <- apply(cur_df[ , cur_samples], 1, function(x) {sum(x == 2, na.rm = T)});
    cur_df$n_total <- cur_df$n_homAnc + cur_df$n_het + cur_df$n_homDer;
    #
    # to save file space, we can calculate this later when we read in the file, instead of printing it now (also has long decimals, waste of file space);
    #cur_df$n_ind <- cur_df$n_total;
    #cur_df$freq_ancestral <- round((2*cur_df$n_homAnc + cur_df$n_het) / (2*cur_df$n_total), 4);
    #cur_df$freq_derived <- round((cur_df$n_het + 2*cur_df$n_homDer) / (2*cur_df$n_total), 4);
    #
    write.table(cur_df[ , c("chr", "pos", "N", "S", "n_homAnc", "n_het", "n_homDer")], paste0("tmp/", cur_pop, "_", allSplits[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }
  
  # combine into a single file;
  # do not print the first line, as it is a header in all files - print the header first as its own separate line;
  # then, in the second command, append using ">>";
  
  system(paste0("cat tmp/", cur_pop, "_aa.txt | head -1 > derivedAlleleFreqs_byPop/derivedAlleleFreqs_", cur_pop, ".txt"));
  system(paste0("tail -n +2 -q tmp/", cur_pop, "* >> derivedAlleleFreqs_byPop/derivedAlleleFreqs_", cur_pop, ".txt"));
  
  # delete the temporary files;
  system(paste0("rm tmp/", cur_pop, "*"));
  
}

#####
#####
#####

##### REPEAT FOR EACH INDIVIDUAL #####

# minor modifications to facilitate a larger number of individuals (use batches for parallelizing) and save file space by not printing chr/pos every time;

##### Process one split-file at a time #####

allGenoFiles <- list.files("genotypes_splitFile/", pattern = "geno_");
allSplits <- gsub("geno_", "", allGenoFiles);

#

# set up batches of 48 to be run at a time;

n_batches <- ceiling(length(allInd) / 48);
for (i in 1:n_batches) {
  if (i == 1) {
    batches = list((48*(i-1) + 1):min(length(allInd), 48 * i));
  } else {
    batches <- append(batches, list((48*(i-1) + 1):min(length(allInd), 48 * i)));
  }
}

#

registerDoParallel(cores=length(allPops));

for (h in 1:length(batches)) {
  
  foreach (k=batches[[h]]) %dopar% {
    
    cur_ind <- allInd[k];

    print(paste0("Working on sample ", cur_ind));
    
    for (i in 1:length(allSplits)) {
      
      print(paste0("Working on split file ", i, " of ", length(allSplits)));
      
      cur_anc <- read.table(paste0("ancestralDerived_splitFile/anc_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
      colnames(cur_anc) <- c("chr", "pos", "ref", "alt", "N", "S", "Macrpy", "Saccja", "Lamidi", "n_outgroup_alleles", "unique_outgroup_allele", "ancestral", "derived");
      
      cur_geno <- read.table(paste0("genotypes_splitFile/geno_", allSplits[i]), header = F, sep = "\t", stringsAsFactors = F);
      colnames(cur_geno) <- c("chr", "pos", "ref", "alt", samples$individual);
      
      # subset to only the genotypes for individuals in the current individual;
      cur_geno <- cur_geno[ , c("chr", "pos", "ref", "alt", cur_ind)];
      
      cur_df <- merge(cur_anc[ , c("chr", "pos", "ref", "alt", "N", "S", "ancestral", "derived")], cur_geno, by = c("chr", "pos", "ref", "alt"), sort = F);
      
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
      # unnecessary to calculate n_homAnc, n_het, n_homDer, n_total as there is only one individual;
      #
      #
      # we only need to output the genotype (genotype has been reformatted to be 0,1,2 format for ancestral/derived), to save file space;
      write.table(cur_df[ , cur_ind], paste0("tmp/", cur_ind, "_", allSplits[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
      
    }
    
    # combine into a single file;
    # do not print the first line, as it is a header in all files - print the header first as its own separate line;
    # then, in the second command, append using ">>";
    
    system(paste0("cat tmp/", cur_ind, "_aa.txt | head -1 > derivedAlleleFreqs_byInd/derivedAlleleFreqs_", cur_ind, ".txt"));
    system(paste0("tail -n +2 -q tmp/", cur_ind, "* >> derivedAlleleFreqs_byInd/derivedAlleleFreqs_", cur_ind, ".txt"));
    
    # delete the temporary files;
    system(paste0("rm tmp/", cur_ind, "*"));
    
  }
  
}
