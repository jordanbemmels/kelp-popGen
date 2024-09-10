##### Script is to calculate realized genetic load in simulated crosses between individuals from all pairwise combinations of populations, including within the same population. This script processes the GERP analysis, and divides the sites into evolutionarily labile vs. evolutionarily conserved sites. #####

# Author use only: GitHub version derived from version 5 of KL14 GERP script hybridF1GTs_v5_methodD_ancestralGERP_KL14.R, 2024/09/10;

#####

### Prerequisites:

# 1) Genotypes of all individuals in each population, pre-printed in [0, 1, 2] indicating the number of copies of the derived allele, previously generated with the getGTs_ancDer_GERP.R. script:

	# genotypes_GERP_byPop/genotypes_derAnc_GERP_[...].txt; # results for each population [...];

# 2) List of individuals that are unrelated to one another, because we do not want simulate crosses between closely-related individuals;
# here, unrelated means that up to first-degree relatives are excluded;
# one individual per line;

	# sampleList_noFirstDegree.txt;

# 3) List of filtered sites to use for calculating realized load, previously calculated using the purging_drift_GERP.R script:

	# global_filteredSitesList_GERP.txt;

##############################
##############################

system("mkdir hybridF1GTs_100reps");

##############################
##############################

# define all populations - these should only be those with a minimum of 3 individuals per population;

allPops <- list.files("genotypes_GERP_byPop");
allPops <- gsub("genotypes_derAnc_GERP_", "", allPops);
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

glbl <- read.table("global_filteredSitesList_GERP.txt", header = T, sep = "\t", stringsAsFactors = F);

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
     
    cur_geno1 <- read.table(paste0("genotypes_GERP_byPop/genotypes_derAnc_GERP_", cur_pop1, ".txt"), sep = "\t", header = T, stringsAsFactors = F);
     colnames(cur_geno1) <- gsub("\\.", "-", colnames(cur_geno1));
      
     cur_geno2 <- read.table(paste0("genotypes_GERP_byPop/genotypes_derAnc_GERP_", cur_pop2, ".txt"), sep = "\t", header = T, stringsAsFactors = F);
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
        
    # split into different dataframes for resampling based on different categories of N and S (GERP results);
    indexes_cur_geno1_LABILE <- which(cur_geno1$N > 0.5 & cur_geno1$S <= 0.5);
    indexes_cur_geno1_CONSERVED <- which(cur_geno1$N > 0.5 & cur_geno1$S > 0.5);
    cur_geno1_LABILE <- cur_geno1[indexes_cur_geno1_LABILE, ];
    cur_geno1_CONSERVED <- cur_geno1[indexes_cur_geno1_CONSERVED, ];

    indexes_cur_geno2_LABILE <- which(cur_geno2$N > 0.5 & cur_geno2$S <= 0.5);
    indexes_cur_geno2_CONSERVED <- which(cur_geno2$N > 0.5 & cur_geno2$S > 0.5);
    cur_geno2_LABILE <- cur_geno2[indexes_cur_geno2_LABILE, ];
    cur_geno2_CONSERVED <- cur_geno2[indexes_cur_geno2_CONSERVED, ];    

    #####
    #####
    
    # LABILE #
    
    # define a set of 100 resampling replicates, picking which individual should be resampled each time;
    # if it is an interpopulation cross there are no restrictions on which individual to use;
    # if it is an intrapopulation cross, we do not want to make selfed individuals, so do not use the same individual twice;
    
    rs100_LABILE <- data.frame(rep = 1:100, pop1 = cur_pop1, pop2 = cur_pop2, parent1 = NA, parent2 = NA, n_loci = NA, rs_nHomAnc = NA, rs_nHet = NA, rs_nHomDer = NA);
    
    for (j in 1:nrow(rs100_LABILE)) {
      if (cur_pop1 == cur_pop2) {
        rs100_LABILE[j, c("parent1", "parent2")] <- sample(colnames(cur_geno1_LABILE)[5:ncol(cur_geno1_LABILE)], size = 2, replace = F);
      } else {
        rs100_LABILE[j , "parent1"] <- sample(colnames(cur_geno1_LABILE)[5:ncol(cur_geno1_LABILE)], size = 1, replace = F);
        rs100_LABILE[j , "parent2"] <- sample(colnames(cur_geno2_LABILE)[5:ncol(cur_geno2_LABILE)], size = 1, replace = F);
      }
    }
    
    #

    # perform the resampling;
    
    for (j in 1:nrow(rs100_LABILE)) {
      #
      # we can only resample if both parents are non-NA genotypes;
      gt_bothParents_LABILE <- cbind(cur_geno1_LABILE[ , rs100_LABILE$parent1[j]], cur_geno2_LABILE[ , rs100_LABILE$parent2[j]]);
      gt_bothParents_LABILE <- gt_bothParents_LABILE[!(is.na(gt_bothParents_LABILE[ , 1])) & !(is.na(gt_bothParents_LABILE[ , 2])), ];
      #
      # ready to resample;
      # note that if genotypes are in [0, 1, 2] format, then the probability of sampling one ancestral allele (0) is (2-x)/2 and the probability of sampling one derived allele (1) is x/2, where x is the genotype;
      alleles_parent1_LABILE <- sapply(gt_bothParents_LABILE[ , 1], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      alleles_parent2_LABILE <- sapply(gt_bothParents_LABILE[ , 2], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      #
      # combine the resampled alleles for both parents into a diploid genotype;
      cur_resampled_GTs_LABILE <- alleles_parent1_LABILE + alleles_parent2_LABILE;
      #
      # count genotypes for output;
      rs100_LABILE$n_loci[j] <- length(cur_resampled_GTs_LABILE);
      rs100_LABILE$rs_nHomAnc[j] <- sum(cur_resampled_GTs_LABILE == 0);
      rs100_LABILE$rs_nHet[j] <- sum(cur_resampled_GTs_LABILE == 1);
      rs100_LABILE$rs_nHomDer[j] <- sum(cur_resampled_GTs_LABILE == 2);
    }

    #####
    #####
    
    # CONSERVED #
    
    # define a set of 100 resampling replicates, picking which individual should be resampled each time;
    # if it is an interpopulation cross there are no restrictions on which individual to use;
    # if it is an intrapopulation cross, we do not want to make selfed individuals, so do not use the same individual twice;
    
    rs100_CONSERVED <- data.frame(rep = 1:100, pop1 = cur_pop1, pop2 = cur_pop2, parent1 = NA, parent2 = NA, n_loci = NA, rs_nHomAnc = NA, rs_nHet = NA, rs_nHomDer = NA);
    
    for (j in 1:nrow(rs100_CONSERVED)) {
      if (cur_pop1 == cur_pop2) {
        rs100_CONSERVED[j, c("parent1", "parent2")] <- sample(colnames(cur_geno1_CONSERVED)[5:ncol(cur_geno1_CONSERVED)], size = 2, replace = F);
      } else {
        rs100_CONSERVED[j , "parent1"] <- sample(colnames(cur_geno1_CONSERVED)[5:ncol(cur_geno1_CONSERVED)], size = 1, replace = F);
        rs100_CONSERVED[j , "parent2"] <- sample(colnames(cur_geno2_CONSERVED)[5:ncol(cur_geno2_CONSERVED)], size = 1, replace = F);
      }
    }
    
    #
    
    # perform the resampling;
    
    for (j in 1:nrow(rs100_CONSERVED)) {
      #
      # we can only resample if both parents are non-NA genotypes;
      gt_bothParents_CONSERVED <- cbind(cur_geno1_CONSERVED[ , rs100_CONSERVED$parent1[j]], cur_geno2_CONSERVED[ , rs100_CONSERVED$parent2[j]]);
      gt_bothParents_CONSERVED <- gt_bothParents_CONSERVED[!(is.na(gt_bothParents_CONSERVED[ , 1])) & !(is.na(gt_bothParents_CONSERVED[ , 2])), ];
      #
      # ready to resample;
      # note that if genotypes are in [0, 1, 2] format, then the probability of sampling one ancestral allele (0) is (2-x)/2 and the probability of sampling one derived allele (1) is x/2, where x is the genotype;
      alleles_parent1_CONSERVED <- sapply(gt_bothParents_CONSERVED[ , 1], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      alleles_parent2_CONSERVED <- sapply(gt_bothParents_CONSERVED[ , 2], function(x) {sample(c(0, 1), size = 1, prob = c((2-x)/2, x/2))});
      #
      # combine the resampled alleles for both parents into a diploid genotype;
      cur_resampled_GTs_CONSERVED <- alleles_parent1_CONSERVED + alleles_parent2_CONSERVED;
      #
      # count genotypes for output;
      rs100_CONSERVED$n_loci[j] <- length(cur_resampled_GTs_CONSERVED);
      rs100_CONSERVED$rs_nHomAnc[j] <- sum(cur_resampled_GTs_CONSERVED == 0);
      rs100_CONSERVED$rs_nHet[j] <- sum(cur_resampled_GTs_CONSERVED == 1);
      rs100_CONSERVED$rs_nHomDer[j] <- sum(cur_resampled_GTs_CONSERVED == 2);
    }
    
    #
    
    write.table(rs100_LABILE, paste0("hybridF1GTs_100reps/rs100_LABILE_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    write.table(rs100_CONSERVED, paste0("hybridF1GTs_100reps/rs100_CONSERVED_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
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
  
  file_rs100_LABILE <- paste0("hybridF1GTs_100reps/rs100_LABILE_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt");
  file_rs100_CONSERVED <- paste0("hybridF1GTs_100reps/rs100_CONSERVED_", all_combos$pop1[i], "_", all_combos$pop2[i], ".txt");
  
  if (file.exists(file_rs100_LABILE)) {
    results_rs100_LABILE <- read.table(file_rs100_LABILE, sep = "\t", header = T, stringsAsFactors = F);
    results_rs100_CONSERVED <- read.table(file_rs100_CONSERVED, sep = "\t", header = T, stringsAsFactors = F);
    #
    cur_df_summary_LABILE <- summarize_resampling(results_rs100_LABILE);
    cur_df_summary_CONSERVED <- summarize_resampling(results_rs100_CONSERVED);
    #
    if (!exists("df_summary_LABILE")) {
      df_summary_LABILE <- cur_df_summary_LABILE;
      df_summary_CONSERVED <- cur_df_summary_CONSERVED;
    } else {
      df_summary_LABILE <- rbind(df_summary_LABILE, cur_df_summary_LABILE);
      df_summary_CONSERVED <- rbind(df_summary_CONSERVED, cur_df_summary_CONSERVED);
    }
  }  
  
}

if(exists("df_summary_LABILE")) {
  write.table(df_summary_LABILE, "hybridF1GTs_100reps_LABILE_summary.txt", sep = "\t", row.names = F, col.names = T, quote = F);  
  write.table(df_summary_CONSERVED, "hybridF1GTs_100reps_CONSERVED_summary.txt", sep = "\t", row.names = F, col.names = T, quote = F);  
}
