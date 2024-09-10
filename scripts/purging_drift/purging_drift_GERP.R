##### Script is to calculate derived allele frequencies overall and by population, to facilitate predictions of purging and genetic drift. This script processes the GERP analysis, and divides the sites into evolutionarily labile vs. evolutionarily conserved sites. #####

# Author use only: GitHub version derived from Part I only (population-level analysis) from version 5 of KL14 GERP script analyzeOutput_GERP_KL14_ancestralGERP_v5.R, 2024/09/10;

#####

### Prerequisites:

# 1) Derived allele frequencies globally across all populations, and separately for each population, both calculated previously using the GERP script get_derivedAlleleFreqs.R:

	# derivedAlleleFreqs.txt; # globally;
	# derivedAlleleFreqs_byPop/derivedAlleleFreqs_[...].txt; # results for each population [...];

#########################
#########################

library(foreach);
library(doParallel);

##############################
##############################

##### LOAD GENOTYPES AND CALCULATE ALLELE FREQUENCIES ACROSS ALL POPULATIONS #####

frq <- read.table("derivedAlleleFreqs.txt", header = T, sep = "\t", stringsAsFactors = F);

frq$n_total <- frq$n_homAnc + frq$n_het + frq$n_homDer;
frq$freq_ancestral <- (2*frq$n_homAnc + frq$n_het) / (2*frq$n_total);
frq$freq_derived <- (frq$n_het + 2*frq$n_homDer) / (2*frq$n_total);
frq$n_derived <- frq$n_het + 2*frq$n_homDer;

##############################
##############################

##### DEFINE SITES OF INTEREST TO BE INCLUDED IN ANALYSES #####

# The number of high-quality sites that should be used in analysis is much lower than the total number of sites printed;
# We should filter these and obtain a final list of sites to analyze in all subsequent analyses, whether globally or per population or per individual;
# Per population or individual may have the sites further drop out due to missing data of course, but we will here define the maximum set of sites to include in downstream analylses;

# Our sites must meet the following criteria:
# 1) ancestral/derived allele must be defined;
# 2) derived allele freq  < 0.5;
# 3) the allele must be variable;
# 4) N (the expected number of mutations on the phylogeny for GERP) must be greater than 0.5, because only if there at least 0.5 expected mutations does it suggest there were enough aligning outgroups in the phylogeny for the GERP score to be meaningful;
# 5) there must be a minimum number of individuals with data - globally, due to our variant calling we likely don't need to implement anything here, but this is worth considering when analyzing individual populations in subsequent steps;

# EXTRA SET:
# 6) we will create two sets of SNPs, one has only the criteria listed above, and a second set (min2) requires a minimum of 2 alleles globally to be included - note that this is implemented in the code here, but the min2 version was not ultimately used in the published manuscript;

# glbl short for global;
glbl <- frq;
nrow(glbl); # 4,172,458;

# filter 1 - ancestral/derived allele must be defined;
glbl <- glbl[!(is.na(glbl$freq_derived)), ];
nrow(glbl); # 3,898,553;

# filter 2 - derived allele must be < 0.5 frequency;
glbl <- glbl[glbl$freq_derived < 0.5, ];
nrow(glbl); # 3,541,802;

# filter 3 - allele must be variable;
glbl <- glbl[glbl$freq_derived != 0 & glbl$freq_derived != 1, ];
nrow(glbl); # 3,467,934;

# filter 4 - N must be greater than 0.5;
glbl <- glbl[glbl$N > 0.5 & !(is.na(glbl$N)), ];
nrow(glbl); # 18,905;

# filter 5 - N must have a minimum  number of individuals;
hist(glbl$n_total); 
min(glbl$n_total); # there is nothing to remove, minimum of 343 individuals;
nrow(glbl); # 18,905;

# filter 6 - for the second set, require a minimum of 2 derived alleles globally;
glbl_min2 <- glbl[glbl$n_derived >= 2, ];
nrow(glbl_min2); # 7,885;
  
# these are our final lists of sites to use in analyses - save to be able to reload later;

write.table(glbl, "global_filteredSitesList_GERP.txt", sep = "\t", row.names = F, col.names = T, quote = F);

write.table(glbl_min2, "global_filteredSitesList_min2_GERP.txt", sep = "\t", row.names = F, col.names = T, quote = F);

##############################
##############################

##### PROCESS BY POPULATION #####

library(foreach);
library(doParallel);

glbl <- read.table("global_filteredSitesList_GERP.txt", header = T, sep = "\t", stringsAsFactors = F);
glbl_min2 <- read.table("global_filteredSitesList_min2_GERP.txt", header = T, sep = "\t", stringsAsFactors = F);

# define all populations;

allPopsFiles <- list.files("derivedAlleleFreqs_byPop");
allPops <- gsub("derivedAlleleFreqs_", "", allPopsFiles);
allPops <- gsub(".txt", "", allPops);

#

# set up batches of 8 to be run at a time;

n_batches <- ceiling(length(allPops) / 8);
for (i in 1:n_batches) {
  if (i == 1) {
    batches = list((8*(i-1) + 1):min(length(allPops), 8 * i));
  } else {
    batches <- append(batches, list((8*(i-1) + 1):min(length(allPops), 8 * i)));
  }
}

# run in parallel;

registerDoParallel(cores=8);

for (k in 1:length(batches)) {
  
  print(paste0("Working on batch ", k, " of ", length(batches)));
  
  foreach (i=batches[[k]]) %dopar% {
    
    cur_pop <- allPops[i];
    cur_pop_file <- allPopsFiles[i];
    
    print(paste0("Working on population ", i, " of ", length(allPops)));
    
    cur_df <- read.table(paste0("derivedAlleleFreqs_byPop/", cur_pop_file), sep = "\t", header = T, stringsAsFactors = F);
    
    # subset to only the filtered sites of interest, for each set;
    # the paste0() makes the sorting with %in% much faster than sorting by two separate columns in each data.frame without pasting;
    cur_df_min0 <- cur_df[paste0(cur_df$chr, "_", cur_df$pos) %in% paste0(glbl$chr, "_", glbl$pos), ];
    cur_df_min2 <- cur_df[paste0(cur_df$chr, "_", cur_df$pos) %in% paste0(glbl_min2$chr, "_", glbl_min2$pos), ];
    
    if (nrow(cur_df_min0) != nrow(glbl)) {
      # skip this population if the number of individuals is not ok;
      print(paste0("ERROR: not all necessary filtered SNPs are present for population ", cur_pop));
      break;
    }
    
    if (nrow(cur_df_min2) != nrow(glbl_min2)) {
      # skip this population if the number of individuals is not ok;
      print(paste0("ERROR: not all necessary filtered SNPs are present for population ", cur_pop));
      break;
    }
    
    #
    
    cur_df_min0$n_total <- cur_df_min0$n_homAnc + cur_df_min0$n_het + cur_df_min0$n_homDer;
    cur_df_min0$freq_ancestral <- (2*cur_df_min0$n_homAnc + cur_df_min0$n_het) / (2*cur_df_min0$n_total);
    cur_df_min0$freq_derived <- (cur_df_min0$n_het + 2*cur_df_min0$n_homDer) / (2*cur_df_min0$n_total);
    
    cur_df_min2$n_total <- cur_df_min2$n_homAnc + cur_df_min2$n_het + cur_df_min2$n_homDer;
    cur_df_min2$freq_ancestral <- (2*cur_df_min2$n_homAnc + cur_df_min2$n_het) / (2*cur_df_min2$n_total);
    cur_df_min2$freq_derived <- (cur_df_min2$n_het + 2*cur_df_min2$n_homDer) / (2*cur_df_min2$n_total);
    
    # some statistics to calculate, for different categories of N and S (different deleteriousness estimates):
    # ACROSS ALL ALLELES:
    # DOES NOT NEED RESAMPLING TO n = 3;
    # mean frequency of derived allele - across all alleles;
    #
    # ACROSS ONLY DERIVED ALLELES THAT ARE PRESENT IN A POPULATION;
    # NEEDS RESAMPLING TO n = 3, SEE BELOW;
    # mean frequency of derived allele - across derived alleles only;
    # proportion of derived alleles that are fixed derived;
    # proportion of derived alleles that are heterozygous in all individuals;
    # proportion of derived alleles that are segregating homozygous (neither entirely fixed derived nor entirely fixed heterozygous);
    
	# the porportion that we observe an allele being fixed is highly dependent upon the number of individuals;
    # thus, we should calculate our statistics by resampling to a minimum of 3 individuals;
    # we can resample only 100 instead of 1000 times, because we only have a maximum of 7 individuals, and 7C3 = 35 so there are only a maximum of 35 unique combinations of individuals anyway, so sampling 1000 times would be overkill and just repeating the same combos of individuals over and over;
    
    # first, we should exclude any alleles genotyped at <3 individuals total - regardless whether we do or do not do any resampling, this is a good idea - note that this means we have different numbers of alleles for each populations and will need to keep track of that too;
    cur_df_min0_min3 <- cur_df_min0[cur_df_min0$n_total >=3, ];
    
    cur_df_min2_min3 <- cur_df_min2[cur_df_min2$n_total >=3, ];
    
    # split into different categories of interest;
    cur_df_min0_LABILE <- cur_df_min0_min3[cur_df_min0_min3$N > 0.5 & cur_df_min0_min3$S <= 0.5, ];
    cur_df_min0_CONSERVED <- cur_df_min0_min3[cur_df_min0_min3$N > 0.5 & cur_df_min0_min3$S > 0.5, ];
    
    cur_df_min2_LABILE <- cur_df_min2_min3[cur_df_min2_min3$N > 0.5 & cur_df_min2_min3$S <= 0.5, ];
    cur_df_min2_CONSERVED <- cur_df_min2_min3[cur_df_min2_min3$N > 0.5 & cur_df_min2_min3$S > 0.5, ];
    
    # STATISTICS CALCULATED ON ALL SITES;
    
    print("Calculating statistics across all sites");
    
    out_df_min0 <- data.frame(pop = cur_pop, LABILE_nSitesAll = NA, LABILE_freq_der_all_mean = NA, LABILE_freq_der_all_sd = NA, LABILE_nSitesDer_total = NA, LABILE_nSitesDer_3ind_mean = NA, LABILE_nSitesDer_3ind_sd = NA, LABILE_freq_der_derOnly_mean = NA, LABILE_freq_der_derOnly_sd = NA, LABILE_prop_der_fixed_mean = NA, LABILE_prop_der_fixed_sd = NA, LABILE_prop_der_allHet_mean = NA, LABILE_prop_der_allHet_sd = NA, LABILE_prop_der_segrHom_mean = NA, LABILE_prop_der_segrHom_sd = NA, LABILE_meanSsites_der_mean = NA, LABILE_meanSsites_der_sd = NA, LABILE_sumSalleles_der_mean = NA, LABILE_sumSalleles_der_sd = NA, LABILE_meanSalleles_der_mean = NA, LABILE_meanSalleles_der_sd = NA, CONSERVED_nSitesAll = NA, CONSERVED_freq_der_all_mean = NA, CONSERVED_freq_der_all_sd = NA, CONSERVED_nSitesDer_total = NA, CONSERVED_nSitesDer_3ind_mean = NA, CONSERVED_nSitesDer_3ind_sd = NA, CONSERVED_freq_der_derOnly_mean = NA, CONSERVED_freq_der_derOnly_sd = NA, CONSERVED_prop_der_fixed_mean = NA, CONSERVED_prop_der_fixed_sd = NA, CONSERVED_prop_der_allHet_mean = NA, CONSERVED_prop_der_allHet_sd = NA, CONSERVED_prop_der_segrHom_mean = NA, CONSERVED_prop_der_segrHom_sd = NA, CONSERVED_meanSsites_der_mean = NA, CONSERVED_meanSsites_der_sd = NA, CONSERVED_sumSalleles_der_mean = NA, CONSERVED_sumSalleles_der_sd = NA, CONSERVED_meanSalleles_der_mean = NA, CONSERVED_meanSalleles_der_sd = NA);

    out_df_min2 <- data.frame(pop = cur_pop, LABILE_nSitesAll = NA, LABILE_freq_der_all_mean = NA, LABILE_freq_der_all_sd = NA, LABILE_nSitesDer_total = NA, LABILE_nSitesDer_3ind_mean = NA, LABILE_nSitesDer_3ind_sd = NA, LABILE_freq_der_derOnly_mean = NA, LABILE_freq_der_derOnly_sd = NA, LABILE_prop_der_fixed_mean = NA, LABILE_prop_der_fixed_sd = NA, LABILE_prop_der_allHet_mean = NA, LABILE_prop_der_allHet_sd = NA, LABILE_prop_der_segrHom_mean = NA, LABILE_prop_der_segrHom_sd = NA, LABILE_meanSsites_der_mean = NA, LABILE_meanSsites_der_sd = NA, LABILE_sumSalleles_der_mean = NA, LABILE_sumSalleles_der_sd = NA, LABILE_meanSalleles_der_mean = NA, LABILE_meanSalleles_der_sd = NA, CONSERVED_nSitesAll = NA, CONSERVED_freq_der_all_mean = NA, CONSERVED_freq_der_all_sd = NA, CONSERVED_nSitesDer_total = NA, CONSERVED_nSitesDer_3ind_mean = NA, CONSERVED_nSitesDer_3ind_sd = NA, CONSERVED_freq_der_derOnly_mean = NA, CONSERVED_freq_der_derOnly_sd = NA, CONSERVED_prop_der_fixed_mean = NA, CONSERVED_prop_der_fixed_sd = NA, CONSERVED_prop_der_allHet_mean = NA, CONSERVED_prop_der_allHet_sd = NA, CONSERVED_prop_der_segrHom_mean = NA, CONSERVED_prop_der_segrHom_sd = NA, CONSERVED_meanSsites_der_mean = NA, CONSERVED_meanSsites_der_sd = NA, CONSERVED_sumSalleles_der_mean = NA, CONSERVED_sumSalleles_der_sd = NA, CONSERVED_meanSalleles_der_mean = NA, CONSERVED_meanSalleles_der_sd = NA);
    
    out_df_min0$LABILE_nSitesAll <- nrow(cur_df_min0_LABILE);
    out_df_min0$LABILE_freq_der_all_mean <- mean(cur_df_min0_LABILE$freq_derived);
    out_df_min0$LABILE_freq_der_all_sd <- sd(cur_df_min0_LABILE$freq_derived);
    
    out_df_min0$CONSERVED_nSitesAll <- nrow(cur_df_min0_CONSERVED);
    out_df_min0$CONSERVED_freq_der_all_mean <- mean(cur_df_min0_CONSERVED$freq_derived);
    out_df_min0$CONSERVED_freq_der_all_sd <- sd(cur_df_min0_CONSERVED$freq_derived);
 
    out_df_min2$LABILE_nSitesAll <- nrow(cur_df_min2_LABILE);
    out_df_min2$LABILE_freq_der_all_mean <- mean(cur_df_min2_LABILE$freq_derived);
    out_df_min2$LABILE_freq_der_all_sd <- sd(cur_df_min2_LABILE$freq_derived);
    
    out_df_min2$CONSERVED_nSitesAll <- nrow(cur_df_min2_CONSERVED);
    out_df_min2$CONSERVED_freq_der_all_mean <- mean(cur_df_min2_CONSERVED$freq_derived);
    out_df_min2$CONSERVED_freq_der_all_sd <- sd(cur_df_min2_CONSERVED$freq_derived);
    
    # STATISTICS CALCULATED ON DERIVED ALLELES ONLY;
    
    print("Calculating statistics across derived alleles only");
    
    # reduce to only derived alleles;
    cur_df_min0_LABILE_derOnly <- cur_df_min0_LABILE[cur_df_min0_LABILE$freq_derived != 0, ];
    cur_df_min0_CONSERVED_derOnly <- cur_df_min0_CONSERVED[cur_df_min0_CONSERVED$freq_derived != 0, ];

    cur_df_min2_LABILE_derOnly <- cur_df_min2_LABILE[cur_df_min2_LABILE$freq_derived != 0, ];
    cur_df_min2_CONSERVED_derOnly <- cur_df_min2_CONSERVED[cur_df_min2_CONSERVED$freq_derived != 0, ];
    
    # calculate statistics that don't require resampling for derived alleles;
    out_df_min0$LABILE_nSitesDer_total <- nrow(cur_df_min0_LABILE_derOnly);
    out_df_min0$CONSERVED_nSitesDer_total <- nrow(cur_df_min0_CONSERVED_derOnly);
 
    out_df_min2$LABILE_nSitesDer_total <- nrow(cur_df_min2_LABILE_derOnly);
    out_df_min2$CONSERVED_nSitesDer_total <- nrow(cur_df_min2_CONSERVED_derOnly);
       
    # set up a dataframe for the results of resampling the derived alleles 100 times;
    
    rs_100_min0 <- data.frame(rep = 1:100, LABILE_nSitesDer_3ind = NA, LABILE_freq_der_derOnly = NA, LABILE_prop_der_fixed = NA, LABILE_prop_der_allHet = NA, LABILE_meanSsites_der = NA, LABILE_sumSalleles_der = NA, LABILE_prop_der_segrHom = NA, LABILE_meanSalleles_der = NA, CONSERVED_nSitesDer_3ind = NA, CONSERVED_freq_der_derOnly = NA, CONSERVED_prop_der_fixed = NA, CONSERVED_prop_der_allHet = NA, CONSERVED_prop_der_segrHom = NA, CONSERVED_meanSsites_der = NA, CONSERVED_sumSalleles_der = NA, CONSERVED_meanSalleles_der = NA);

    rs_100_min2 <- data.frame(rep = 1:100, LABILE_nSitesDer_3ind = NA, LABILE_freq_der_derOnly = NA, LABILE_prop_der_fixed = NA, LABILE_prop_der_allHet = NA, LABILE_meanSsites_der = NA, LABILE_sumSalleles_der = NA, LABILE_prop_der_segrHom = NA, LABILE_meanSalleles_der = NA, CONSERVED_nSitesDer_3ind = NA, CONSERVED_freq_der_derOnly = NA, CONSERVED_prop_der_fixed = NA, CONSERVED_prop_der_allHet = NA, CONSERVED_prop_der_segrHom = NA, CONSERVED_meanSsites_der = NA, CONSERVED_sumSalleles_der = NA, CONSERVED_meanSalleles_der = NA);
    
    # resample derived alleles only and calculate statistics;
    
    # z must be supplied as the columns c("n_homAnc", "n_het", "n_homDer", "S") - we resample only the first three;
    resample_by_line <- function(z) {
      z_rs <- sample(rep(x = c(0, 1, 2), times = z[1:3]), size = 3, replace = F);
      z_rs <- c(z_rs, z[4]); # add the "S" vector position to the resampled results;
      z_rs;
    }
        
    #
    
    if (nrow(cur_df_min0_LABILE_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: min0, LABILE");
      
      for (j in 1:100) {
        
        rs_LABILE <- as.data.frame(t(apply(cur_df_min0_LABILE_derOnly[ , c("n_homAnc", "n_het", "n_homDer", "S")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_LABILE$freq_der <- (rs_LABILE$V1 + rs_LABILE$V2 + rs_LABILE$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_LABILE <- rs_LABILE[rs_LABILE$freq_der > 0, ];
        
        rs_LABILE$is_fixed_der <- rs_LABILE$freq_der == 1;
        rs_LABILE$is_allHet <- apply(rs_LABILE[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_LABILE$is_segrHom <- (rs_LABILE$is_fixed_der == F) & (rs_LABILE$is_allHet == F);

        rs_LABILE$Ssites_der <- rs_LABILE$S;
        rs_LABILE$Salleles_der <- apply(rs_LABILE[ , 1:4], 1, function(x) {sum(x[1:3]) * x[4]}); # each allele copy (0, 1, 2 per individual) is used to calculate the total number of derived alleles, then multiplied by S, to get the total "S" load contribution summed across all alleles;
                
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100_min0$LABILE_nSitesDer_3ind[j] <- nrow(rs_LABILE);
        rs_100_min0$LABILE_freq_der_derOnly[j] <- mean(rs_LABILE$freq_der);
        rs_100_min0$LABILE_prop_der_fixed[j] <- mean(rs_LABILE$is_fixed_der);
        rs_100_min0$LABILE_prop_der_allHet[j] <- mean(rs_LABILE$is_allHet);
        rs_100_min0$LABILE_prop_der_segrHom[j] <- mean(rs_LABILE$is_segrHom);
        rs_100_min0$LABILE_meanSsites_der[j] <- mean(rs_LABILE$Ssites_der);
        rs_100_min0$LABILE_sumSalleles_der[j] <- sum(rs_LABILE$Salleles_der);
        rs_100_min0$LABILE_meanSalleles_der[j] <- mean(rs_LABILE$Salleles_der);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df_min0$LABILE_nSitesDer_3ind_mean <- mean(rs_100_min0$LABILE_nSitesDer_3ind);
      out_df_min0$LABILE_nSitesDer_3ind_sd <- sd(rs_100_min0$LABILE_nSitesDer_3ind);
      out_df_min0$LABILE_freq_der_derOnly_mean <- mean(rs_100_min0$LABILE_freq_der_derOnly);
      out_df_min0$LABILE_freq_der_derOnly_sd <- sd(rs_100_min0$LABILE_freq_der_derOnly);
      out_df_min0$LABILE_prop_der_fixed_mean <- mean(rs_100_min0$LABILE_prop_der_fixed);
      out_df_min0$LABILE_prop_der_fixed_sd <- sd(rs_100_min0$LABILE_prop_der_fixed);
      out_df_min0$LABILE_prop_der_allHet_mean <- mean(rs_100_min0$LABILE_prop_der_allHet);
      out_df_min0$LABILE_prop_der_allHet_sd <- sd(rs_100_min0$LABILE_prop_der_allHet);
      out_df_min0$LABILE_prop_der_segrHom_mean <- mean(rs_100_min0$LABILE_prop_der_segrHom);
      out_df_min0$LABILE_prop_der_segrHom_sd <- sd(rs_100_min0$LABILE_prop_der_segrHom);
      out_df_min0$LABILE_meanSsites_der_mean <- mean(rs_100_min0$LABILE_meanSsites_der);
      out_df_min0$LABILE_meanSsites_der_sd <- sd(rs_100_min0$LABILE_meanSsites_der);
      out_df_min0$LABILE_sumSalleles_der_mean <- mean(rs_100_min0$LABILE_sumSalleles_der);
      out_df_min0$LABILE_sumSalleles_der_sd <- sd(rs_100_min0$LABILE_sumSalleles_der);
      out_df_min0$LABILE_meanSalleles_der_mean <- mean(rs_100_min0$LABILE_meanSalleles_der);
      out_df_min0$LABILE_meanSalleles_der_sd <- sd(rs_100_min0$LABILE_meanSalleles_der);
      out_df_min0$LABILE_meanSalleles_der_mean <- mean(rs_100_min0$LABILE_meanSalleles_der);
      out_df_min0$LABILE_meanSalleles_der_sd <- sd(rs_100_min0$LABILE_meanSalleles_der);
      
    } 
    
    #
    
    if (nrow(cur_df_min0_CONSERVED_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: min0, CONSERVED");
      
      for (j in 1:100) {
        
        rs_CONSERVED <- as.data.frame(t(apply(cur_df_min0_CONSERVED_derOnly[ , c("n_homAnc", "n_het", "n_homDer", "S")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_CONSERVED$freq_der <- (rs_CONSERVED$V1 + rs_CONSERVED$V2 + rs_CONSERVED$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_CONSERVED <- rs_CONSERVED[rs_CONSERVED$freq_der > 0, ];
        
        rs_CONSERVED$is_fixed_der <- rs_CONSERVED$freq_der == 1;
        rs_CONSERVED$is_allHet <- apply(rs_CONSERVED[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_CONSERVED$is_segrHom <- (rs_CONSERVED$is_fixed_der == F) & (rs_CONSERVED$is_allHet == F);
        
        rs_CONSERVED$Ssites_der <- rs_CONSERVED$S;
        rs_CONSERVED$Salleles_der <- apply(rs_CONSERVED[ , 1:4], 1, function(x) {sum(x[1:3]) * x[4]}); # each allele copy (0, 1, 2 per individual) is used to calculate the total number of derived alleles, then multiplied by S, to get the total "S" load contribution summed across all alleles;
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100_min0$CONSERVED_nSitesDer_3ind[j] <- nrow(rs_CONSERVED);
        rs_100_min0$CONSERVED_freq_der_derOnly[j] <- mean(rs_CONSERVED$freq_der);
        rs_100_min0$CONSERVED_prop_der_fixed[j] <- mean(rs_CONSERVED$is_fixed_der);
        rs_100_min0$CONSERVED_prop_der_allHet[j] <- mean(rs_CONSERVED$is_allHet);
        rs_100_min0$CONSERVED_prop_der_segrHom[j] <- mean(rs_CONSERVED$is_segrHom);
        rs_100_min0$CONSERVED_meanSsites_der[j] <- mean(rs_CONSERVED$Ssites_der);
        rs_100_min0$CONSERVED_sumSalleles_der[j] <- sum(rs_CONSERVED$Salleles_der);
        rs_100_min0$CONSERVED_meanSalleles_der[j] <- mean(rs_CONSERVED$Salleles_der);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df_min0$CONSERVED_nSitesDer_3ind_mean <- mean(rs_100_min0$CONSERVED_nSitesDer_3ind);
      out_df_min0$CONSERVED_nSitesDer_3ind_sd <- sd(rs_100_min0$CONSERVED_nSitesDer_3ind);
      out_df_min0$CONSERVED_freq_der_derOnly_mean <- mean(rs_100_min0$CONSERVED_freq_der_derOnly);
      out_df_min0$CONSERVED_freq_der_derOnly_sd <- sd(rs_100_min0$CONSERVED_freq_der_derOnly);
      out_df_min0$CONSERVED_prop_der_fixed_mean <- mean(rs_100_min0$CONSERVED_prop_der_fixed);
      out_df_min0$CONSERVED_prop_der_fixed_sd <- sd(rs_100_min0$CONSERVED_prop_der_fixed);
      out_df_min0$CONSERVED_prop_der_allHet_mean <- mean(rs_100_min0$CONSERVED_prop_der_allHet);
      out_df_min0$CONSERVED_prop_der_allHet_sd <- sd(rs_100_min0$CONSERVED_prop_der_allHet);
      out_df_min0$CONSERVED_prop_der_segrHom_mean <- mean(rs_100_min0$CONSERVED_prop_der_segrHom);
      out_df_min0$CONSERVED_prop_der_segrHom_sd <- sd(rs_100_min0$CONSERVED_prop_der_segrHom);
      out_df_min0$CONSERVED_meanSsites_der_mean <- mean(rs_100_min0$CONSERVED_meanSsites_der);
      out_df_min0$CONSERVED_meanSsites_der_sd <- sd(rs_100_min0$CONSERVED_meanSsites_der);
      out_df_min0$CONSERVED_sumSalleles_der_mean <- mean(rs_100_min0$CONSERVED_sumSalleles_der);
      out_df_min0$CONSERVED_sumSalleles_der_sd <- sd(rs_100_min0$CONSERVED_sumSalleles_der);
      out_df_min0$CONSERVED_meanSalleles_der_mean <- mean(rs_100_min0$CONSERVED_meanSalleles_der);
      out_df_min0$CONSERVED_meanSalleles_der_sd <- sd(rs_100_min0$CONSERVED_meanSalleles_der);
      out_df_min0$CONSERVED_meanSalleles_der_mean <- mean(rs_100_min0$CONSERVED_meanSalleles_der);
      out_df_min0$CONSERVED_meanSalleles_der_sd <- sd(rs_100_min0$CONSERVED_meanSalleles_der);
      
    } 
    
    #
    
    if (nrow(cur_df_min2_LABILE_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: min2, LABILE");
      
      for (j in 1:100) {
        
        rs_LABILE <- as.data.frame(t(apply(cur_df_min2_LABILE_derOnly[ , c("n_homAnc", "n_het", "n_homDer", "S")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_LABILE$freq_der <- (rs_LABILE$V1 + rs_LABILE$V2 + rs_LABILE$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_LABILE <- rs_LABILE[rs_LABILE$freq_der > 0, ];
        
        rs_LABILE$is_fixed_der <- rs_LABILE$freq_der == 1;
        rs_LABILE$is_allHet <- apply(rs_LABILE[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_LABILE$is_segrHom <- (rs_LABILE$is_fixed_der == F) & (rs_LABILE$is_allHet == F);
        
        rs_LABILE$Ssites_der <- rs_LABILE$S;
        rs_LABILE$Salleles_der <- apply(rs_LABILE[ , 1:4], 1, function(x) {sum(x[1:3]) * x[4]}); # each allele copy (0, 1, 2 per individual) is used to calculate the total number of derived alleles, then multiplied by S, to get the total "S" load contribution summed across all alleles;
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100_min2$LABILE_nSitesDer_3ind[j] <- nrow(rs_LABILE);
        rs_100_min2$LABILE_freq_der_derOnly[j] <- mean(rs_LABILE$freq_der);
        rs_100_min2$LABILE_prop_der_fixed[j] <- mean(rs_LABILE$is_fixed_der);
        rs_100_min2$LABILE_prop_der_allHet[j] <- mean(rs_LABILE$is_allHet);
        rs_100_min2$LABILE_prop_der_segrHom[j] <- mean(rs_LABILE$is_segrHom);
        rs_100_min2$LABILE_meanSsites_der[j] <- mean(rs_LABILE$Ssites_der);
        rs_100_min2$LABILE_sumSalleles_der[j] <- sum(rs_LABILE$Salleles_der);
        rs_100_min2$LABILE_meanSalleles_der[j] <- mean(rs_LABILE$Salleles_der);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df_min2$LABILE_nSitesDer_3ind_mean <- mean(rs_100_min2$LABILE_nSitesDer_3ind);
      out_df_min2$LABILE_nSitesDer_3ind_sd <- sd(rs_100_min2$LABILE_nSitesDer_3ind);
      out_df_min2$LABILE_freq_der_derOnly_mean <- mean(rs_100_min2$LABILE_freq_der_derOnly);
      out_df_min2$LABILE_freq_der_derOnly_sd <- sd(rs_100_min2$LABILE_freq_der_derOnly);
      out_df_min2$LABILE_prop_der_fixed_mean <- mean(rs_100_min2$LABILE_prop_der_fixed);
      out_df_min2$LABILE_prop_der_fixed_sd <- sd(rs_100_min2$LABILE_prop_der_fixed);
      out_df_min2$LABILE_prop_der_allHet_mean <- mean(rs_100_min2$LABILE_prop_der_allHet);
      out_df_min2$LABILE_prop_der_allHet_sd <- sd(rs_100_min2$LABILE_prop_der_allHet);
      out_df_min2$LABILE_prop_der_segrHom_mean <- mean(rs_100_min2$LABILE_prop_der_segrHom);
      out_df_min2$LABILE_prop_der_segrHom_sd <- sd(rs_100_min2$LABILE_prop_der_segrHom);
      out_df_min2$LABILE_meanSsites_der_mean <- mean(rs_100_min2$LABILE_meanSsites_der);
      out_df_min2$LABILE_meanSsites_der_sd <- sd(rs_100_min2$LABILE_meanSsites_der);
      out_df_min2$LABILE_sumSalleles_der_mean <- mean(rs_100_min2$LABILE_sumSalleles_der);
      out_df_min2$LABILE_sumSalleles_der_sd <- sd(rs_100_min2$LABILE_sumSalleles_der);
      out_df_min2$LABILE_meanSalleles_der_mean <- mean(rs_100_min2$LABILE_meanSalleles_der);
      out_df_min2$LABILE_meanSalleles_der_sd <- sd(rs_100_min2$LABILE_meanSalleles_der);
      out_df_min2$LABILE_meanSalleles_der_mean <- mean(rs_100_min2$LABILE_meanSalleles_der);
      out_df_min2$LABILE_meanSalleles_der_sd <- sd(rs_100_min2$LABILE_meanSalleles_der);
      
    } 
    
    #
    
    if (nrow(cur_df_min2_CONSERVED_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: min2, CONSERVED");
      
      for (j in 1:100) {
        
        rs_CONSERVED <- as.data.frame(t(apply(cur_df_min2_CONSERVED_derOnly[ , c("n_homAnc", "n_het", "n_homDer", "S")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_CONSERVED$freq_der <- (rs_CONSERVED$V1 + rs_CONSERVED$V2 + rs_CONSERVED$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_CONSERVED <- rs_CONSERVED[rs_CONSERVED$freq_der > 0, ];
        
        rs_CONSERVED$is_fixed_der <- rs_CONSERVED$freq_der == 1;
        rs_CONSERVED$is_allHet <- apply(rs_CONSERVED[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_CONSERVED$is_segrHom <- (rs_CONSERVED$is_fixed_der == F) & (rs_CONSERVED$is_allHet == F);
        
        rs_CONSERVED$Ssites_der <- rs_CONSERVED$S;
        rs_CONSERVED$Salleles_der <- apply(rs_CONSERVED[ , 1:4], 1, function(x) {sum(x[1:3]) * x[4]}); # each allele copy (0, 1, 2 per individual) is used to calculate the total number of derived alleles, then multiplied by S, to get the total "S" load contribution summed across all alleles;
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100_min2$CONSERVED_nSitesDer_3ind[j] <- nrow(rs_CONSERVED);
        rs_100_min2$CONSERVED_freq_der_derOnly[j] <- mean(rs_CONSERVED$freq_der);
        rs_100_min2$CONSERVED_prop_der_fixed[j] <- mean(rs_CONSERVED$is_fixed_der);
        rs_100_min2$CONSERVED_prop_der_allHet[j] <- mean(rs_CONSERVED$is_allHet);
        rs_100_min2$CONSERVED_prop_der_segrHom[j] <- mean(rs_CONSERVED$is_segrHom);
        rs_100_min2$CONSERVED_meanSsites_der[j] <- mean(rs_CONSERVED$Ssites_der);
        rs_100_min2$CONSERVED_sumSalleles_der[j] <- sum(rs_CONSERVED$Salleles_der);
        rs_100_min2$CONSERVED_meanSalleles_der[j] <- mean(rs_CONSERVED$Salleles_der);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df_min2$CONSERVED_nSitesDer_3ind_mean <- mean(rs_100_min2$CONSERVED_nSitesDer_3ind);
      out_df_min2$CONSERVED_nSitesDer_3ind_sd <- sd(rs_100_min2$CONSERVED_nSitesDer_3ind);
      out_df_min2$CONSERVED_freq_der_derOnly_mean <- mean(rs_100_min2$CONSERVED_freq_der_derOnly);
      out_df_min2$CONSERVED_freq_der_derOnly_sd <- sd(rs_100_min2$CONSERVED_freq_der_derOnly);
      out_df_min2$CONSERVED_prop_der_fixed_mean <- mean(rs_100_min2$CONSERVED_prop_der_fixed);
      out_df_min2$CONSERVED_prop_der_fixed_sd <- sd(rs_100_min2$CONSERVED_prop_der_fixed);
      out_df_min2$CONSERVED_prop_der_allHet_mean <- mean(rs_100_min2$CONSERVED_prop_der_allHet);
      out_df_min2$CONSERVED_prop_der_allHet_sd <- sd(rs_100_min2$CONSERVED_prop_der_allHet);
      out_df_min2$CONSERVED_prop_der_segrHom_mean <- mean(rs_100_min2$CONSERVED_prop_der_segrHom);
      out_df_min2$CONSERVED_prop_der_segrHom_sd <- sd(rs_100_min2$CONSERVED_prop_der_segrHom);
      out_df_min2$CONSERVED_meanSsites_der_mean <- mean(rs_100_min2$CONSERVED_meanSsites_der);
      out_df_min2$CONSERVED_meanSsites_der_sd <- sd(rs_100_min2$CONSERVED_meanSsites_der);
      out_df_min2$CONSERVED_sumSalleles_der_mean <- mean(rs_100_min2$CONSERVED_sumSalleles_der);
      out_df_min2$CONSERVED_sumSalleles_der_sd <- sd(rs_100_min2$CONSERVED_sumSalleles_der);
      out_df_min2$CONSERVED_meanSalleles_der_mean <- mean(rs_100_min2$CONSERVED_meanSalleles_der);
      out_df_min2$CONSERVED_meanSalleles_der_sd <- sd(rs_100_min2$CONSERVED_meanSalleles_der);
      out_df_min2$CONSERVED_meanSalleles_der_mean <- mean(rs_100_min2$CONSERVED_meanSalleles_der);
      out_df_min2$CONSERVED_meanSalleles_der_sd <- sd(rs_100_min2$CONSERVED_meanSalleles_der);
      
    } 
    
    #
 
    # write output;
    
    write.table(out_df_min0, paste0("outStats_GERP/outStats_min0_", cur_pop, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    write.table(out_df_min2, paste0("outStats_GERP/outStats_min2_", cur_pop, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }
  
}

#

# combine outfiles into a single file;
# do not print the first line, as it is a header in all files - print the header first as its own separate line;
# then, in the second command, append using ">>";

system(paste0("cat outStats_GERP/outStats_min0_", allPops[1], ".txt | head -1 > outStats_GERP/outStats_allPops_min0.txt"));
system(paste0("tail -n +2 -q outStats_GERP/outStats_min0_*.txt >> outStats_GERP/outStats_allPops_min0.txt"));

system(paste0("cat outStats_GERP/outStats_min2_", allPops[1], ".txt | head -1 > outStats_GERP/outStats_allPops_min2.txt"));
system(paste0("tail -n +2 -q outStats_GERP/outStats_min2_*.txt >> outStats_GERP/outStats_allPops_min2.txt"));
