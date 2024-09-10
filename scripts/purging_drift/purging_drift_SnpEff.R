##### Script is to calculate derived allele frequencies overall and by population, to facilitate predictions of purging and genetic drift. This script processes the SnpEff analysis, and divides the sites into modifier and low-, medium-, and high-impact sites. #####

# Author use only: GitHub version derived from Part I only (population-level analysis) from version Anc1 of KL14 SnpEff script analyzeOutput_snpEff_KL14_KRI_ancestralGERP_vAnc1.R, 2024/09/10;

#####

### Prerequisites:

# 1) Derived allele frequencies globally across all populations, and separately for each population, both calculated previously using the SnpEff script getDerivedAlleleFreqs_snpEff.R:

	# derivedAlleleFreqs_snpEff_allPops.txt.txt; # globally;
	# derivedAlleleFreqs_snpEff_byPop/derivedAlleleFreqs_snpEff_[...].txt; # results for each population [...];

#########################
#########################

library(foreach);
library(doParallel);

##############################
##############################

##### INSPECT RESULTS ACROSS ALL POPULATIONS - PRELIMINARY #####

frq <- read.table("derivedAlleleFreqs_snpEff_allPops.txt", header = T, sep = "\t", stringsAsFactors = F);

frq$n_total <- frq$n_homAnc + frq$n_het + frq$n_homDer;
frq$freq_ancestral <- (2*frq$n_homAnc + frq$n_het) / (2*frq$n_total);
frq$freq_derived <- (frq$n_het + 2*frq$n_homDer) / (2*frq$n_total);

##############################
##############################

##### DEFINE SITES OF INTEREST TO BE INCLUDED IN ANALYSES #####

# The number of high-quality sites that should be used in analysis is much lower than the total number of sites printed;
# We should filter these and obtain a final list of sites to analyze in all subsequent analyses, whether globally or per population or per individual;
# Per population or individual may have the sites further drop out due to missing data of course, but we will here define the maximum set of sites to include in downstream analylses;

# Our sites must meet the following criteria:
# 1) ancestral/derived allele must be defined - this was already taken care of for SnpEff (but wasn't for GERP at this point);
# 2) derived allele freq  < 0.5;
# 3) the allele must be variable;
# 4) [criterion about GERP scores, not relevant to SnpEff];
# 5) there must be a minimum number of individuals with data - globally, due to our variant calling we likely don't need to implement anything here, but this is worth considering when analyzing individual populations in subsequent steps;

# glbl short for global;
glbl <- frq;
nrow(glbl); # 3,820,266;

# filter 1 - ancestral/derived allele must be defined;
glbl <- glbl[!(is.na(glbl$freq_derived)), ];
nrow(glbl); # 3,820,266;

# filter 2 - derived allele must be < 0.5 frequency;
glbl <- glbl[glbl$freq_derived < 0.5, ];
nrow(glbl); # 3,467,934;

# filter 3 - allele must be variable;
glbl <- glbl[glbl$freq_derived != 0 & glbl$freq_derived != 1, ];
nrow(glbl); # 3,467,934;

# filter 4 - N must be greater than zero;
# NA;

# filter 5 - N must have a minimum  number of individuals;
hist(glbl$n_total); 
min(glbl$n_total); # there is nothing to remove, minimum of 341 individuals;
nrow(glbl); # 3,467,934;
  
# this is our final list of sites to use in analyses - save to be able to reload later;

write.table(glbl, "global_filteredSitesList_snpEff.txt", sep = "\t", row.names = F, col.names = T, quote = F);

##############################
##############################

##### PROCESS BY POPULATION #####

library(foreach);
library(doParallel);

glbl <- read.table("global_filteredSitesList_snpEff.txt", header = T, sep = "\t", stringsAsFactors = F);

# define all populations;

allPopsFiles <- list.files("derivedAlleleFreqs_snpEff_byPop");
allPops <- gsub("derivedAlleleFreqs_snpEff_", "", allPopsFiles);
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
    
    cur_df <- read.table(paste0("derivedAlleleFreqs_snpEff_byPop/", cur_pop_file), sep = "\t", header = T, stringsAsFactors = F);
    
    # subset to only the filtered sites of interest;
    # the paste0() makes the sorting with %in% much faster than sorting by two separate columns in each data.frame without pasting;
    cur_df <- cur_df[paste0(cur_df$chr, "_", cur_df$pos) %in% paste0(glbl$chr, "_", glbl$pos), ];
    
    if (nrow(cur_df) != nrow(glbl)) {
      # skip this population if the number of individuals is not ok;
      print(paste0("ERROR: not all necessary filtered SNPs are present for population ", cur_pop));
      break;
    }
    
    #
    
    cur_df$n_total <- cur_df$n_homAnc + cur_df$n_het + cur_df$n_homDer;
    cur_df$freq_ancestral <- (2*cur_df$n_homAnc + cur_df$n_het) / (2*cur_df$n_total);
    cur_df$freq_derived <- (cur_df$n_het + 2*cur_df$n_homDer) / (2*cur_df$n_total);
    
    # some statistics to calculate, for different categories putative_impact (different deleteriousness estimates):
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
    cur_df_min3 <- cur_df[cur_df$n_total >=3, ];
    
    # split into different categories of interest;
    cur_df_MODIFIER <- cur_df_min3[cur_df_min3$putative_impact == "MODIFIER", ];
    cur_df_LOW <- cur_df_min3[cur_df_min3$putative_impact == "LOW", ];
    cur_df_MODERATE <- cur_df_min3[cur_df_min3$putative_impact == "MODERATE", ];
    cur_df_HIGH <- cur_df_min3[cur_df_min3$putative_impact == "HIGH", ];
    
    # STATISTICS CALCULATED ON MODIFIER SITES;
    
    print("Calculating statistics across MODIFIER sites");
    
    out_df <- data.frame(pop = cur_pop, MODIFIER_nSitesAll = NA, MODIFIER_freq_der_all_mean = NA, MODIFIER_freq_der_all_sd = NA, MODIFIER_nSitesDer_total = NA, MODIFIER_nSitesDer_3ind_mean = NA, MODIFIER_nSitesDer_3ind_sd = NA, MODIFIER_freq_der_derOnly_mean = NA, MODIFIER_freq_der_derOnly_sd = NA, MODIFIER_prop_der_fixed_mean = NA, MODIFIER_prop_der_fixed_sd = NA, MODIFIER_prop_der_allHet_mean = NA, MODIFIER_prop_der_allHet_sd = NA, MODIFIER_prop_der_segrHom_mean = NA, MODIFIER_prop_der_segrHom_sd = NA, LOW_nSitesAll = NA, LOW_freq_der_all_mean = NA, LOW_freq_der_all_sd = NA, LOW_nSitesDer_total = NA, LOW_nSitesDer_3ind_mean = NA, LOW_nSitesDer_3ind_sd = NA, LOW_freq_der_derOnly_mean = NA, LOW_freq_der_derOnly_sd = NA, LOW_prop_der_fixed_mean = NA, LOW_prop_der_fixed_sd = NA, LOW_prop_der_allHet_mean = NA, LOW_prop_der_allHet_sd = NA, LOW_prop_der_segrHom_mean = NA, LOW_prop_der_segrHom_sd = NA, MODERATE_nSitesAll = NA, MODERATE_freq_der_all_mean = NA, MODERATE_freq_der_all_sd = NA, MODERATE_nSitesDer_total = NA, MODERATE_nSitesDer_3ind_mean = NA, MODERATE_nSitesDer_3ind_sd = NA, MODERATE_freq_der_derOnly_mean = NA, MODERATE_freq_der_derOnly_sd = NA, MODERATE_prop_der_fixed_mean = NA, MODERATE_prop_der_fixed_sd = NA, MODERATE_prop_der_allHet_mean = NA, MODERATE_prop_der_allHet_sd = NA, MODERATE_prop_der_segrHom_mean = NA, MODERATE_prop_der_segrHom_sd = NA, HIGH_nSitesAll = NA, HIGH_freq_der_all_mean = NA, HIGH_freq_der_all_sd = NA, HIGH_nSitesDer_total = NA, HIGH_nSitesDer_3ind_mean = NA, HIGH_nSitesDer_3ind_sd = NA, HIGH_freq_der_derOnly_mean = NA, HIGH_freq_der_derOnly_sd = NA, HIGH_prop_der_fixed_mean = NA, HIGH_prop_der_fixed_sd = NA, HIGH_prop_der_allHet_mean = NA, HIGH_prop_der_allHet_sd = NA, HIGH_prop_der_segrHom_mean = NA, HIGH_prop_der_segrHom_sd = NA);
    
    out_df$MODIFIER_nSitesAll <- nrow(cur_df_MODIFIER);
    out_df$MODIFIER_freq_der_all_mean <- mean(cur_df_MODIFIER$freq_derived);
    out_df$MODIFIER_freq_der_all_sd <- sd(cur_df_MODIFIER$freq_derived);
    
    out_df$LOW_nSitesAll <- nrow(cur_df_LOW);
    out_df$LOW_freq_der_all_mean <- mean(cur_df_LOW$freq_derived);
    out_df$LOW_freq_der_all_sd <- sd(cur_df_LOW$freq_derived);
    
    out_df$MODERATE_nSitesAll <- nrow(cur_df_MODERATE);
    out_df$MODERATE_freq_der_all_mean <- mean(cur_df_MODERATE$freq_derived);
    out_df$MODERATE_freq_der_all_sd <- sd(cur_df_MODERATE$freq_derived);
    
    out_df$HIGH_nSitesAll <- nrow(cur_df_HIGH);
    out_df$HIGH_freq_der_all_mean <- mean(cur_df_HIGH$freq_derived);
    out_df$HIGH_freq_der_all_sd <- sd(cur_df_HIGH$freq_derived);
    
    # STATISTICS CALCULATED ON DERIVED ALLELES ONLY;
    
    print("Calculating statistics across derived alleles only");
    
    # reduce to only derived alleles;
    cur_df_MODIFIER_derOnly <- cur_df_MODIFIER[cur_df_MODIFIER$freq_derived != 0, ];
    cur_df_LOW_derOnly <- cur_df_LOW[cur_df_LOW$freq_derived != 0, ];
    cur_df_MODERATE_derOnly <- cur_df_MODERATE[cur_df_MODERATE$freq_derived != 0, ];
    cur_df_HIGH_derOnly <- cur_df_HIGH[cur_df_HIGH$freq_derived != 0, ];
    
    # calculate statistics that don't require resampling for derived alleles;
    out_df$MODIFIER_nSitesDer_total <- nrow(cur_df_MODIFIER_derOnly);
    out_df$LOW_nSitesDer_total <- nrow(cur_df_LOW_derOnly);
    out_df$MODERATE_nSitesDer_total <- nrow(cur_df_MODERATE_derOnly);
    out_df$HIGH_nSitesDer_total <- nrow(cur_df_HIGH_derOnly);
    
    # set up a dataframe for the results of resampling the derived alleles 100 times;
    
    rs_100 <- data.frame(rep = 1:100, MODIFIER_nSitesDer_3ind = NA, MODIFIER_freq_der_derOnly = NA, MODIFIER_prop_der_fixed = NA, MODIFIER_prop_der_allHet = NA, MODIFIER_prop_der_segrHom = NA, MODIFIER_meanSsites_der = NA, MODIFIER_sumSalleles_der = NA, MODIFIER_meanSalleles_der = NA, LOW_nSitesDer_3ind = NA, LOW_freq_der_derOnly = NA, LOW_prop_der_fixed = NA, LOW_prop_der_allHet = NA, LOW_prop_der_segrHom = NA, LOW_meanSsites_der = NA, LOW_sumSalleles_der = NA, LOW_meanSalleles_der = NA, MODERATE_nSitesDer_3ind = NA, MODERATE_freq_der_derOnly = NA, MODERATE_prop_der_fixed = NA, MODERATE_prop_der_allHet = NA, MODERATE_meanSsites_der = NA, MODERATE_sumSalleles_der = NA, MODERATE_prop_der_segrHom = NA, MODERATE_meanSalleles_der = NA, HIGH_nSitesDer_3ind = NA, HIGH_freq_der_derOnly = NA, HIGH_prop_der_fixed = NA, HIGH_prop_der_allHet = NA, HIGH_prop_der_segrHom = NA, HIGH_meanSsites_der = NA, HIGH_sumSalleles_der = NA, HIGH_meanSalleles_der = NA);

    # resample derived alleles only and calculate statistics;
    
    # z must be supplied as the columns c("n_homAnc", "n_het", "n_homDer") - we resample only the first three;
    resample_by_line <- function(z) {
      z_rs <- sample(rep(x = c(0, 1, 2), times = z[1:3]), size = 3, replace = F);
      z_rs;
    }
    
    if (nrow(cur_df_MODIFIER_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: MODIFIER");
      
      for (j in 1:100) {
        
        rs_MODIFIER <- as.data.frame(t(apply(cur_df_MODIFIER_derOnly[ , c("n_homAnc", "n_het", "n_homDer")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_MODIFIER$freq_der <- (rs_MODIFIER$V1 + rs_MODIFIER$V2 + rs_MODIFIER$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_MODIFIER <- rs_MODIFIER[rs_MODIFIER$freq_der > 0, ];
        
        rs_MODIFIER$is_fixed_der <- rs_MODIFIER$freq_der == 1;
        rs_MODIFIER$is_allHet <- apply(rs_MODIFIER[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_MODIFIER$is_segrHom <- (rs_MODIFIER$is_fixed_der == F) & (rs_MODIFIER$is_allHet == F);
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100$MODIFIER_nSitesDer_3ind[j] <- nrow(rs_MODIFIER);
        rs_100$MODIFIER_freq_der_derOnly[j] <- mean(rs_MODIFIER$freq_der);
        rs_100$MODIFIER_prop_der_fixed[j] <- mean(rs_MODIFIER$is_fixed_der);
        rs_100$MODIFIER_prop_der_allHet[j] <- mean(rs_MODIFIER$is_allHet);
        rs_100$MODIFIER_prop_der_segrHom[j] <- mean(rs_MODIFIER$is_segrHom);

      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df$MODIFIER_nSitesDer_3ind_mean <- mean(rs_100$MODIFIER_nSitesDer_3ind);
      out_df$MODIFIER_nSitesDer_3ind_sd <- sd(rs_100$MODIFIER_nSitesDer_3ind);
      out_df$MODIFIER_freq_der_derOnly_mean <- mean(rs_100$MODIFIER_freq_der_derOnly);
      out_df$MODIFIER_freq_der_derOnly_sd <- sd(rs_100$MODIFIER_freq_der_derOnly);
      out_df$MODIFIER_prop_der_fixed_mean <- mean(rs_100$MODIFIER_prop_der_fixed);
      out_df$MODIFIER_prop_der_fixed_sd <- sd(rs_100$MODIFIER_prop_der_fixed);
      out_df$MODIFIER_prop_der_allHet_mean <- mean(rs_100$MODIFIER_prop_der_allHet);
      out_df$MODIFIER_prop_der_allHet_sd <- sd(rs_100$MODIFIER_prop_der_allHet);
      out_df$MODIFIER_prop_der_segrHom_mean <- mean(rs_100$MODIFIER_prop_der_segrHom);
      out_df$MODIFIER_prop_der_segrHom_sd <- sd(rs_100$MODIFIER_prop_der_segrHom);

    }
    
    #
    
    if (nrow(cur_df_LOW_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: LOW");
      
      for (j in 1:100) {
        
        rs_LOW <- as.data.frame(t(apply(cur_df_LOW_derOnly[ , c("n_homAnc", "n_het", "n_homDer")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_LOW$freq_der <- (rs_LOW$V1 + rs_LOW$V2 + rs_LOW$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_LOW <- rs_LOW[rs_LOW$freq_der > 0, ];
        
        rs_LOW$is_fixed_der <- rs_LOW$freq_der == 1;
        rs_LOW$is_allHet <- apply(rs_LOW[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_LOW$is_segrHom <- (rs_LOW$is_fixed_der == F) & (rs_LOW$is_allHet == F);
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100$LOW_nSitesDer_3ind[j] <- nrow(rs_LOW);
        rs_100$LOW_freq_der_derOnly[j] <- mean(rs_LOW$freq_der);
        rs_100$LOW_prop_der_fixed[j] <- mean(rs_LOW$is_fixed_der);
        rs_100$LOW_prop_der_allHet[j] <- mean(rs_LOW$is_allHet);
        rs_100$LOW_prop_der_segrHom[j] <- mean(rs_LOW$is_segrHom);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df$LOW_nSitesDer_3ind_mean <- mean(rs_100$LOW_nSitesDer_3ind);
      out_df$LOW_nSitesDer_3ind_sd <- sd(rs_100$LOW_nSitesDer_3ind);
      out_df$LOW_freq_der_derOnly_mean <- mean(rs_100$LOW_freq_der_derOnly);
      out_df$LOW_freq_der_derOnly_sd <- sd(rs_100$LOW_freq_der_derOnly);
      out_df$LOW_prop_der_fixed_mean <- mean(rs_100$LOW_prop_der_fixed);
      out_df$LOW_prop_der_fixed_sd <- sd(rs_100$LOW_prop_der_fixed);
      out_df$LOW_prop_der_allHet_mean <- mean(rs_100$LOW_prop_der_allHet);
      out_df$LOW_prop_der_allHet_sd <- sd(rs_100$LOW_prop_der_allHet);
      out_df$LOW_prop_der_segrHom_mean <- mean(rs_100$LOW_prop_der_segrHom);
      out_df$LOW_prop_der_segrHom_sd <- sd(rs_100$LOW_prop_der_segrHom);
      
    }
    
    #
    
    if (nrow(cur_df_MODERATE_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: MODERATE");
      
      for (j in 1:100) {
        
        rs_MODERATE <- as.data.frame(t(apply(cur_df_MODERATE_derOnly[ , c("n_homAnc", "n_het", "n_homDer")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_MODERATE$freq_der <- (rs_MODERATE$V1 + rs_MODERATE$V2 + rs_MODERATE$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_MODERATE <- rs_MODERATE[rs_MODERATE$freq_der > 0, ];
        
        rs_MODERATE$is_fixed_der <- rs_MODERATE$freq_der == 1;
        rs_MODERATE$is_allHet <- apply(rs_MODERATE[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_MODERATE$is_segrHom <- (rs_MODERATE$is_fixed_der == F) & (rs_MODERATE$is_allHet == F);
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100$MODERATE_nSitesDer_3ind[j] <- nrow(rs_MODERATE);
        rs_100$MODERATE_freq_der_derOnly[j] <- mean(rs_MODERATE$freq_der);
        rs_100$MODERATE_prop_der_fixed[j] <- mean(rs_MODERATE$is_fixed_der);
        rs_100$MODERATE_prop_der_allHet[j] <- mean(rs_MODERATE$is_allHet);
        rs_100$MODERATE_prop_der_segrHom[j] <- mean(rs_MODERATE$is_segrHom);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df$MODERATE_nSitesDer_3ind_mean <- mean(rs_100$MODERATE_nSitesDer_3ind);
      out_df$MODERATE_nSitesDer_3ind_sd <- sd(rs_100$MODERATE_nSitesDer_3ind);
      out_df$MODERATE_freq_der_derOnly_mean <- mean(rs_100$MODERATE_freq_der_derOnly);
      out_df$MODERATE_freq_der_derOnly_sd <- sd(rs_100$MODERATE_freq_der_derOnly);
      out_df$MODERATE_prop_der_fixed_mean <- mean(rs_100$MODERATE_prop_der_fixed);
      out_df$MODERATE_prop_der_fixed_sd <- sd(rs_100$MODERATE_prop_der_fixed);
      out_df$MODERATE_prop_der_allHet_mean <- mean(rs_100$MODERATE_prop_der_allHet);
      out_df$MODERATE_prop_der_allHet_sd <- sd(rs_100$MODERATE_prop_der_allHet);
      out_df$MODERATE_prop_der_segrHom_mean <- mean(rs_100$MODERATE_prop_der_segrHom);
      out_df$MODERATE_prop_der_segrHom_sd <- sd(rs_100$MODERATE_prop_der_segrHom);
      
    }
    
    #
    
    if (nrow(cur_df_HIGH_derOnly) > 0) {
      
      print("Performing resampling 100 times for derived alleles: HIGH");
      
      for (j in 1:100) {
        
        rs_HIGH <- as.data.frame(t(apply(cur_df_HIGH_derOnly[ , c("n_homAnc", "n_het", "n_homDer")], 1, resample_by_line)));
        
        # we need to calculate the derived allele frequency and EXCLUDE any sites that no longer have any derived allele present;
        rs_HIGH$freq_der <- (rs_HIGH$V1 + rs_HIGH$V2 + rs_HIGH$V3) / 6; # we can calculate quickly this way because 0 = 0 copies of derived allele, 1 = 1 copy, 2 = 2 copies, and there are always 6 alleles total;
        rs_HIGH <- rs_HIGH[rs_HIGH$freq_der > 0, ];
        
        rs_HIGH$is_fixed_der <- rs_HIGH$freq_der == 1;
        rs_HIGH$is_allHet <- apply(rs_HIGH[ , 1:3], 1, function(x) {sum(x == 1) == 3});
        rs_HIGH$is_segrHom <- (rs_HIGH$is_fixed_der == F) & (rs_HIGH$is_allHet == F);
        
        # ready to calculate the values for traits of interest for this resampling iteration;
        rs_100$HIGH_nSitesDer_3ind[j] <- nrow(rs_HIGH);
        rs_100$HIGH_freq_der_derOnly[j] <- mean(rs_HIGH$freq_der);
        rs_100$HIGH_prop_der_fixed[j] <- mean(rs_HIGH$is_fixed_der);
        rs_100$HIGH_prop_der_allHet[j] <- mean(rs_HIGH$is_allHet);
        rs_100$HIGH_prop_der_segrHom[j] <- mean(rs_HIGH$is_segrHom);
        
      }
      
      # calculate the values for the traits that should be derived from our set of 100 sampling replicates;
      # we are placing this here within the {} structure so that we only update the values if the resampling was done, and if it wasn't done that means there weren't any sites left to resample and we can leave the values as NA;
      out_df$HIGH_nSitesDer_3ind_mean <- mean(rs_100$HIGH_nSitesDer_3ind);
      out_df$HIGH_nSitesDer_3ind_sd <- sd(rs_100$HIGH_nSitesDer_3ind);
      out_df$HIGH_freq_der_derOnly_mean <- mean(rs_100$HIGH_freq_der_derOnly);
      out_df$HIGH_freq_der_derOnly_sd <- sd(rs_100$HIGH_freq_der_derOnly);
      out_df$HIGH_prop_der_fixed_mean <- mean(rs_100$HIGH_prop_der_fixed);
      out_df$HIGH_prop_der_fixed_sd <- sd(rs_100$HIGH_prop_der_fixed);
      out_df$HIGH_prop_der_allHet_mean <- mean(rs_100$HIGH_prop_der_allHet);
      out_df$HIGH_prop_der_allHet_sd <- sd(rs_100$HIGH_prop_der_allHet);
      out_df$HIGH_prop_der_segrHom_mean <- mean(rs_100$HIGH_prop_der_segrHom);
      out_df$HIGH_prop_der_segrHom_sd <- sd(rs_100$HIGH_prop_der_segrHom);
      
    }
    
    # write output;
    
    write.table(out_df, paste0("outStats_snpEff/outStats_snpEff_", cur_pop, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }
  
}

#

# combine outfiles into a single file;
# do not print the first line, as it is a header in all files - print the header first as its own separate line;
# then, in the second command, append using ">>";

system(paste0("cat outStats_snpEff/outStats_snpEff_", allPops[1], ".txt | head -1 > outStats_snpEff/outStats_allPops_snpEff.txt"));
system(paste0("tail -n +2 -q outStats_snpEff/outStats_snpEff_*.txt >> outStats_snpEff/outStats_allPops_snpEff.txt"));
