##### Script is to calculate realized genetic load in simulated crosses between individuals from all pairwise combinations of populations, including within the same population. This script processes the SnpEff analysis, and divides the sites into modifier and low-, moderate-, and high-impact sites. #####

# Author use only: GitHub version derived from version 5 of KL14 SnpEff script hybridF1GTs_vAnc1_methodD_ancestralGERP_KL14_KRI, 2024/09/10;

#####

### Prerequisites:

# 1) Genotypes of all individuals in each population, pre-printed in [0, 1, 2] indicating the number of copies of the derived allele, previously generated with the getGTs_ancDer_SnpEff.R. script:

	# genotypes_snpEff_byPop/genotypes_derAnc_snpEff_[...].txt; # results for each population [...];

# 2) List of individuals that are unrelated to one another, because we do not want simulate crosses between closely-related individuals;
# here, unrelated means that up to first-degree relatives are excluded;
# one individual per line;

	# sampleList_noFirstDegree.txt;

# 3) List of filtered sites to use for calculating realized load, previously calculated using the purging_drift_SnpEff.R script:

	# global_filteredSitesList_snpEff.txt;

##############################
##############################

system("mkdir hybridF1GTs_100reps");

##############################
##############################

# define all populations - these should only be those with a minimum of 3 individuals per population;

allPops <- list.files("genotypes_snpEff_byPop");
allPops <- gsub("genotypes_derAnc_snpEff_", "", allPops);
allPops <- gsub(".txt", "", allPops);

# set up all combinations of populations;

all_combos <- data.frame(pop1 = NA, pop2 = NA);
for (i in 1:length(allPops)) {
  for (j in 1:length(allPops)) {
    if (i <= j) {
      new_pop1 <- allPops[i];
      new_pop2 <- allPops[j];
      all_combos <- rbind(all_combos, c(new_pop1, new_pop2));
    }
  }
}

all_combos <- all_combos[!(is.na(all_combos$pop1)), ];

##############################
##############################

# get a list of individuals with first-degree relatives excluded - for this analysis we do not want to include close relatives because we want the realized load to be for matings between unrelated individuals only within each population, otherwise the statistics may be biased to partly reflect the effects of inbreeding rather than differential fixation of alleles by genetic drift;
# we will use this list below to ensure that we only sample individuals on this list;

noFirstDegree <- read.table("sampleList_noFirstDegree.txt", sep = "\t", header = F, stringsAsFactors = F);

# user may need to edit this line if the individual is indicated differently relative to the sample name;
noFirstDegree$V1 <- gsub("KL..\\.", "", noFirstDegree$V1);

##############################
##############################

# load the list of sites to include;

glbl <- read.table("global_filteredSitesList_snpEff.txt", header = T, sep = "\t", stringsAsFactors = F);

##############################
##############################

library(foreach);
library(doParallel);

# set up 48 to be run at a time;

#n_batches <- 48;
length_batches <- 48;
n_batches <- ceiling(nrow(all_combos) / length_batches);
for (i in 1:n_batches) {
  if (i == 1) {
    batches = list((length_batches*(i-1) + 1):min(nrow(all_combos), length_batches * i));
  } else {
    batches <- append(batches, list((length_batches*(i-1) + 1):min(nrow(all_combos), length_batches * i)));
  }
}

#

registerDoParallel(cores=48);

for (k in 1:length(batches)) {
  
  foreach (i=batches[[k]]) %dopar% {

    print(paste0("Working on combination ", i, " of ", nrow(all_combos)));
    
    cur_pop1 <- all_combos$pop1[i];
    cur_pop2 <- all_combos$pop2[i];
    
     # load the genotypes of the current populations;
     
    cur_geno1 <- read.table(paste0("genotypes_snpEff_byPop/genotypes_derAnc_snpEff_", cur_pop1, ".txt"), sep = "\t", header = T, stringsAsFactors = F);
     colnames(cur_geno1) <- gsub("\\.", "-", colnames(cur_geno1));
      
     cur_geno2 <- read.table(paste0("genotypes_snpEff_byPop/genotypes_derAnc_snpEff_", cur_pop2, ".txt"), sep = "\t", header = T, stringsAsFactors = F);
     colnames(cur_geno2) <- gsub("\\.", "-", colnames(cur_geno2));
    
    #
    
    # subset to only the filtered sites of interest;
    # the paste0() makes the sorting with %in% much faster than sorting by two separate columns in each data.frame without pasting;
    cur_geno1 <- cur_geno1[paste0(cur_geno1$chr, cur_geno1$pos) %in% paste0(glbl$chr, glbl$pos), ];
    cur_geno2 <- cur_geno2[paste0(cur_geno2$chr, cur_geno2$pos) %in% paste0(glbl$chr, glbl$pos), ];    
    
    #
    
    # remove individuals not present in the noFirstDegree dataframe;
    cur_geno1 <- cur_geno1[ , c(1:4, which(colnames(cur_geno1) %in% noFirstDegree$V1))];
    cur_geno2 <- cur_geno2[ , c(1:4, which(colnames(cur_geno2) %in% noFirstDegree$V1))];

    # split into different dataframes for resampling based on different snpEff categories;
    indexes_cur_geno1_MODIFIER <- which(cur_geno1$putative_impact == "MODIFIER");
    indexes_cur_geno1_LOW <- which(cur_geno1$putative_impact == "LOW");
    indexes_cur_geno1_MODERATE <- which(cur_geno1$putative_impact == "MODERATE");
    indexes_cur_geno1_HIGH <- which(cur_geno1$putative_impact == "HIGH");
    
    cur_geno1_MODIFIER <- cur_geno1[indexes_cur_geno1_MODIFIER, ];
    cur_geno1_LOW <- cur_geno1[indexes_cur_geno1_LOW, ];
    cur_geno1_MODERATE <- cur_geno1[indexes_cur_geno1_MODERATE, ];
    cur_geno1_HIGH <- cur_geno1[indexes_cur_geno1_HIGH, ];
 
    indexes_cur_geno2_MODIFIER <- which(cur_geno2$putative_impact == "MODIFIER");
    indexes_cur_geno2_LOW <- which(cur_geno2$putative_impact == "LOW");
    indexes_cur_geno2_MODERATE <- which(cur_geno2$putative_impact == "MODERATE");
    indexes_cur_geno2_HIGH <- which(cur_geno2$putative_impact == "HIGH");
    
    cur_geno2_MODIFIER <- cur_geno2[indexes_cur_geno2_MODIFIER, ];
    cur_geno2_LOW <- cur_geno2[indexes_cur_geno2_LOW, ];
    cur_geno2_MODERATE <- cur_geno2[indexes_cur_geno2_MODERATE, ];
    cur_geno2_HIGH <- cur_geno2[indexes_cur_geno2_HIGH, ];
    
    #####
    #####
    
    # DO NOT RUN MODIFIER AS WE LIKELY WON'T USE IT AND IT IS TOO SLOW#
    # MODIFIER #
    
    # define a set of 100 resampling replicates, picking which individual should be resampled each time;
    # if it is an interpopulation cross there are no restrictions on which individual to use;
    # if it is an intrapopulation cross, we do not want to make selfed individuals, so do not use the same individual twice;
    
#    rs100_MODIFIER <- data.frame(rep = 1:100, pop1 = cur_pop1, pop2 = cur_pop2, parent1 = NA, parent2 = NA, n_loci = NA, rs_nHomAnc = NA, rs_nHet = NA, rs_nHomDer = NA);
#    
#    for (j in 1:nrow(rs100_MODIFIER)) {
#      if (cur_pop1 == cur_pop2) {
#        rs100_MODIFIER[j, c("parent1", "parent2")] <- sample(colnames(cur_geno1_MODIFIER)[5:ncol(cur_geno1_MODIFIER)], size = 2, replace = F);
#      } else {
#        rs100_MODIFIER[j , "parent1"] <- sample(colnames(cur_geno1_MODIFIER)[5:ncol(cur_geno1_MODIFIER)], size = 1, replace = F);
#        rs100_MODIFIER[j , "parent2"] <- sample(colnames(cur_geno2_MODIFIER)[5:ncol(cur_geno2_MODIFIER)], size = 1, replace = F);
#      }
#    }
    
    #
    
    # perform the resampling;
    
#    for (j in 1:nrow(rs100_MODIFIER)) {
      #
      # we can only resample if both parents are non-NA genotypes;
#      gt_bothParents_MODIFIER <- cbind(cur_geno1_MODIFIER[ , rs100_MODIFIER$parent1[j]], cur_geno2_MODIFIER[ , rs100_MODIFIER$parent2[j]]);
#      gt_bothParents_MODIFIER <- gt_bothParents_MODIFIER[!(is.na(gt_bothParents_MODIFIER[ , 1])) & !(is.na(gt_bothParents_MODIFIER[ , 2])), ];
      #
      # ready to resample;
      # note that if genotypes are in [0, 1, 2] format, then the probability of sampling one ancestral allele (0) is (2-x)/2 and the probability of sampling one derived allele (1) is x/2, where x is the genotype;
#      alleles_parent1_MODIFIER <- sapply(gt_bothParents_MODIFIER[ , 1], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
#      alleles_parent2_MODIFIER <- sapply(gt_bothParents_MODIFIER[ , 2], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      #
      # combine the resampled alleles for both parents into a diploid genotype;
#      cur_resampled_GTs_MODIFIER <- alleles_parent1_MODIFIER + alleles_parent2_MODIFIER;
      #
      # count genotypes for output;
#      rs100_MODIFIER$n_loci[j] <- length(cur_resampled_GTs_MODIFIER);
#      rs100_MODIFIER$rs_nHomAnc[j] <- sum(cur_resampled_GTs_MODIFIER == 0);
#      rs100_MODIFIER$rs_nHet[j] <- sum(cur_resampled_GTs_MODIFIER == 1);
#      rs100_MODIFIER$rs_nHomDer[j] <- sum(cur_resampled_GTs_MODIFIER == 2);
#    }    
    
    #####
    #####
    
    # LOW #
    
    # define a set of 100 resampling replicates, picking which individual should be resampled each time;
    # if it is an interpopulation cross there are no restrictions on which individual to use;
    # if it is an intrapopulation cross, we do not want to make selfed individuals, so do not use the same individual twice;
    
    rs100_LOW <- data.frame(rep = 1:100, pop1 = cur_pop1, pop2 = cur_pop2, parent1 = NA, parent2 = NA, n_loci = NA, rs_nHomAnc = NA, rs_nHet = NA, rs_nHomDer = NA);
    
    for (j in 1:nrow(rs100_LOW)) {
      if (cur_pop1 == cur_pop2) {
        rs100_LOW[j, c("parent1", "parent2")] <- sample(colnames(cur_geno1_LOW)[5:ncol(cur_geno1_LOW)], size = 2, replace = F);
      } else {
        rs100_LOW[j , "parent1"] <- sample(colnames(cur_geno1_LOW)[5:ncol(cur_geno1_LOW)], size = 1, replace = F);
        rs100_LOW[j , "parent2"] <- sample(colnames(cur_geno2_LOW)[5:ncol(cur_geno2_LOW)], size = 1, replace = F);
      }
    }
    
    #
    
    # perform the resampling;
    
    for (j in 1:nrow(rs100_LOW)) {
      #
      # we can only resample if both parents are non-NA genotypes;
      gt_bothParents_LOW <- cbind(cur_geno1_LOW[ , rs100_LOW$parent1[j]], cur_geno2_LOW[ , rs100_LOW$parent2[j]]);
      gt_bothParents_LOW <- gt_bothParents_LOW[!(is.na(gt_bothParents_LOW[ , 1])) & !(is.na(gt_bothParents_LOW[ , 2])), ];
      #
      # ready to resample;
      # note that if genotypes are in [0, 1, 2] format, then the probability of sampling one ancestral allele (0) is (2-x)/2 and the probability of sampling one derived allele (1) is x/2, where x is the genotype;
      alleles_parent1_LOW <- sapply(gt_bothParents_LOW[ , 1], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      alleles_parent2_LOW <- sapply(gt_bothParents_LOW[ , 2], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      #
      # combine the resampled alleles for both parents into a diploid genotype;
      cur_resampled_GTs_LOW <- alleles_parent1_LOW + alleles_parent2_LOW;
      #
      # count genotypes for output;
      rs100_LOW$n_loci[j] <- length(cur_resampled_GTs_LOW);
      rs100_LOW$rs_nHomAnc[j] <- sum(cur_resampled_GTs_LOW == 0);
      rs100_LOW$rs_nHet[j] <- sum(cur_resampled_GTs_LOW == 1);
      rs100_LOW$rs_nHomDer[j] <- sum(cur_resampled_GTs_LOW == 2);
    }
    
    #####
    #####
    
    # MODERATE #
    
    # define a set of 100 resampling replicates, picking which individual should be resampled each time;
    # if it is an interpopulation cross there are no restrictions on which individual to use;
    # if it is an intrapopulation cross, we do not want to make selfed individuals, so do not use the same individual twice;
    
    rs100_MODERATE <- data.frame(rep = 1:100, pop1 = cur_pop1, pop2 = cur_pop2, parent1 = NA, parent2 = NA, n_loci = NA, rs_nHomAnc = NA, rs_nHet = NA, rs_nHomDer = NA);
    
    for (j in 1:nrow(rs100_MODERATE)) {
      if (cur_pop1 == cur_pop2) {
        rs100_MODERATE[j, c("parent1", "parent2")] <- sample(colnames(cur_geno1_MODERATE)[5:ncol(cur_geno1_MODERATE)], size = 2, replace = F);
      } else {
        rs100_MODERATE[j , "parent1"] <- sample(colnames(cur_geno1_MODERATE)[5:ncol(cur_geno1_MODERATE)], size = 1, replace = F);
        rs100_MODERATE[j , "parent2"] <- sample(colnames(cur_geno2_MODERATE)[5:ncol(cur_geno2_MODERATE)], size = 1, replace = F);
      }
    }
    
    #
    
    # perform the resampling;
    
    for (j in 1:nrow(rs100_MODERATE)) {
      #
      # we can only resample if both parents are non-NA genotypes;
      gt_bothParents_MODERATE <- cbind(cur_geno1_MODERATE[ , rs100_MODERATE$parent1[j]], cur_geno2_MODERATE[ , rs100_MODERATE$parent2[j]]);
      gt_bothParents_MODERATE <- gt_bothParents_MODERATE[!(is.na(gt_bothParents_MODERATE[ , 1])) & !(is.na(gt_bothParents_MODERATE[ , 2])), ];
      #
      # ready to resample;
      # note that if genotypes are in [0, 1, 2] format, then the probability of sampling one ancestral allele (0) is (2-x)/2 and the probability of sampling one derived allele (1) is x/2, where x is the genotype;
      alleles_parent1_MODERATE <- sapply(gt_bothParents_MODERATE[ , 1], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      alleles_parent2_MODERATE <- sapply(gt_bothParents_MODERATE[ , 2], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      #
      # combine the resampled alleles for both parents into a diploid genotype;
      cur_resampled_GTs_MODERATE <- alleles_parent1_MODERATE + alleles_parent2_MODERATE;
      #
      # count genotypes for output;
      rs100_MODERATE$n_loci[j] <- length(cur_resampled_GTs_MODERATE);
      rs100_MODERATE$rs_nHomAnc[j] <- sum(cur_resampled_GTs_MODERATE == 0);
      rs100_MODERATE$rs_nHet[j] <- sum(cur_resampled_GTs_MODERATE == 1);
      rs100_MODERATE$rs_nHomDer[j] <- sum(cur_resampled_GTs_MODERATE == 2);
    }
    
    #####
    #####
    
    # HIGH #
    
    # define a set of 100 resampling replicates, picking which individual should be resampled each time;
    # if it is an interpopulation cross there are no restrictions on which individual to use;
    # if it is an intrapopulation cross, we do not want to make selfed individuals, so do not use the same individual twice;
    
    rs100_HIGH <- data.frame(rep = 1:100, pop1 = cur_pop1, pop2 = cur_pop2, parent1 = NA, parent2 = NA, n_loci = NA, rs_nHomAnc = NA, rs_nHet = NA, rs_nHomDer = NA);
    
    for (j in 1:nrow(rs100_HIGH)) {
      if (cur_pop1 == cur_pop2) {
        rs100_HIGH[j, c("parent1", "parent2")] <- sample(colnames(cur_geno1_HIGH)[5:ncol(cur_geno1_HIGH)], size = 2, replace = F);
      } else {
        rs100_HIGH[j , "parent1"] <- sample(colnames(cur_geno1_HIGH)[5:ncol(cur_geno1_HIGH)], size = 1, replace = F);
        rs100_HIGH[j , "parent2"] <- sample(colnames(cur_geno2_HIGH)[5:ncol(cur_geno2_HIGH)], size = 1, replace = F);
      }
    }
    
    #
    
    # perform the resampling;
    
    for (j in 1:nrow(rs100_HIGH)) {
      #
      # we can only resample if both parents are non-NA genotypes;
      gt_bothParents_HIGH <- cbind(cur_geno1_HIGH[ , rs100_HIGH$parent1[j]], cur_geno2_HIGH[ , rs100_HIGH$parent2[j]]);
      gt_bothParents_HIGH <- gt_bothParents_HIGH[!(is.na(gt_bothParents_HIGH[ , 1])) & !(is.na(gt_bothParents_HIGH[ , 2])), ];
      #
      # ready to resample;
      # note that if genotypes are in [0, 1, 2] format, then the probability of sampling one ancestral allele (0) is (2-x)/2 and the probability of sampling one derived allele (1) is x/2, where x is the genotype;
      alleles_parent1_HIGH <- sapply(gt_bothParents_HIGH[ , 1], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      alleles_parent2_HIGH <- sapply(gt_bothParents_HIGH[ , 2], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      #
      # combine the resampled alleles for both parents into a diploid genotype;
      cur_resampled_GTs_HIGH <- alleles_parent1_HIGH + alleles_parent2_HIGH;
      #
      # count genotypes for output;
      rs100_HIGH$n_loci[j] <- length(cur_resampled_GTs_HIGH);
      rs100_HIGH$rs_nHomAnc[j] <- sum(cur_resampled_GTs_HIGH == 0);
      rs100_HIGH$rs_nHet[j] <- sum(cur_resampled_GTs_HIGH == 1);
      rs100_HIGH$rs_nHomDer[j] <- sum(cur_resampled_GTs_HIGH == 2);
    }
    
    #
    
    #write.table(rs100_MODIFIER, paste0("hybridF1GTs_100reps/rs100_MODIFIER_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    write.table(rs100_LOW, paste0("hybridF1GTs_100reps/rs100_LOW_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    write.table(rs100_MODERATE, paste0("hybridF1GTs_100reps/rs100_MODERATE_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    write.table(rs100_HIGH, paste0("hybridF1GTs_100reps/rs100_HIGH_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }  
  
}

##############################
##############################

# function to summarize across all population combinations;

summarize_resampling <- function(rs_results) {
  
  rs_results$rs_prop_homAnc <- rs_results$rs_nHomAnc / rs_results$n_loci;
  rs_results$rs_prop_het <- rs_results$rs_nHet / rs_results$n_loci;
  rs_results$rs_prop_homDer <- rs_results$rs_nHomDer / rs_results$n_loci;
  
  rs_results$rs_freqAnc <- (2*rs_results$rs_nHomAnc + rs_results$rs_nHet) / (2*rs_results$n_loci);
  rs_results$rs_freqDer <- (rs_results$rs_nHet + 2*rs_results$rs_nHomDer) / (2*rs_results$n_loci);
  
  rs_results_means <- as.data.frame(t(apply(rs_results[ , c("rep", "n_loci", "rs_prop_homAnc", "rs_prop_het", "rs_prop_homDer", "rs_freqAnc", "rs_freqDer")], 2, mean)));
  
  rs_results_means$pop1 <- unique(rs_results$pop1);
  rs_results_means$pop2 <- unique(rs_results$pop2);
  
  return(rs_results_means);
  
}

# summarize across all individuals;

for (i in 1:nrow(all_combos)) {
  
  file_rs100_LOW <- paste0("hybridF1GTs_100reps/rs100_LOW_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt");
  file_rs100_MODERATE <- paste0("hybridF1GTs_100reps/rs100_MODERATE_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt");
  file_rs100_HIGH <- paste0("hybridF1GTs_100reps/rs100_HIGH_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt");
  
  if (file.exists(file_rs100_LOW)) {
    results_rs100_LOW <- read.table(file_rs100_LOW, sep = "\t", header = T, stringsAsFactors = F);
    results_rs100_MODERATE <- read.table(file_rs100_MODERATE, sep = "\t", header = T, stringsAsFactors = F);
    results_rs100_HIGH <- read.table(file_rs100_HIGH, sep = "\t", header = T, stringsAsFactors = F);
    #
    cur_df_summary_LOW <- summarize_resampling(results_rs100_LOW);
    cur_df_summary_MODERATE <- summarize_resampling(results_rs100_MODERATE);
    cur_df_summary_HIGH <- summarize_resampling(results_rs100_HIGH);
        #
    if (!exists("df_summary_LOW")) {
      df_summary_LOW <- cur_df_summary_LOW;
      df_summary_MODERATE <- cur_df_summary_MODERATE;
      df_summary_HIGH <- cur_df_summary_HIGH;
          } else {
      df_summary_LOW <- rbind(df_summary_LOW, cur_df_summary_LOW);
      df_summary_MODERATE <- rbind(df_summary_MODERATE, cur_df_summary_MODERATE);
      df_summary_HIGH <- rbind(df_summary_HIGH, cur_df_summary_HIGH);
          }
  }  
  
}

if(exists("df_summary_LOW")) {
  write.table(df_summary_LOW, "hybridF1GTs_100reps_LOW_summary.txt", sep = "\t", row.names = F, col.names = T, quote = F);  
  write.table(df_summary_MODERATE, "hybridF1GTs_100reps_MODERATE_summary.txt", sep = "\t", row.names = F, col.names = T, quote = F);  
  write.table(df_summary_HIGH, "hybridF1GTs_100reps_HIGH_summary.txt", sep = "\t", row.names = F, col.names = T, quote = F);  

}
