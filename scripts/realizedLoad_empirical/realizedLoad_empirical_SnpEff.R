##### Script is to calculate realized genetic load in each individual. This script processes the SnpEff analysis, and divides the sites into modifier and low-, moderate-, and high-impact sites. #####

# Author use only: GitHub version derived from Part II only (individual-level analysis) from version Anc1 of KL14 SnpEff script analyzeOutput_snpEff_KL14_KRI_ancestralGERP_vAnc1.R, 2024/09/10;

#####

### Prerequisites:

# 1) Derived allele frequencies separately for each population, and separately for each individual, all calculated previously using the SnpEff script getDerivedAlleleFreqs_snpEff.R:

	# derivedAlleleFreqs_snpEff_byPop/derivedAlleleFreqs_snpEff_[...].txt; # results for each population [...];
	# derivedAlleleFreqs_snpEff_byInd/derivedAlleleFreqs_snpEff_[...].txt; # results for each individual [...];

# 2) List of filtered sites to use for calculating realized load, previously calculated using the purging_drift_SnpEff.R script:

	# global_filteredSitesList_snpEff.txt;
	
#########################
#########################

library(foreach);
library(doParallel);

##############################
##############################

##### PROCESS BY INDIVIDUAL #####

library(foreach);
library(doParallel);

# re-load our list of good sites;
glbl <- read.table("global_filteredSitesList_snpEff.txt", header = T, sep = "\t", stringsAsFactors = F);

# re-load our list of all sites, so that we know which line is which site in the individual files - this was not outputted, in order to save space;
# just load the first population as an example, here NL-BA-01 - we don't need all populations;
sitesInOrder <- read.table("derivedAlleleFreqs_snpEff_byPop/derivedAlleleFreqs_snpEff_NL-BA-01.txt", sep = "\t", header = T, stringsAsFactors = F);
sitesInOrder <- sitesInOrder[ , c("chr", "pos")];

goodSitesIndex <- which(paste0(sitesInOrder$chr, "_", sitesInOrder$pos) %in% paste0(glbl$chr, "_", glbl$pos));
length(goodSitesIndex) == nrow(glbl);

# define all individuals;

allIndFiles <- list.files("derivedAlleleFreqs_snpEff_byInd");
allInd <- gsub("derivedAlleleFreqs_snpEff_", "", allIndFiles);
allInd <- gsub(".txt", "", allInd);

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

# run in parallel;

registerDoParallel(cores=48);

for (k in 1:length(batches)) {
  
  print(paste0("Working on batch ", k, " of ", length(batches)));
  
  foreach (i=batches[[k]]) %dopar% {
    
    cur_ind <- allInd[i];
    cur_ind_file <- allIndFiles[i];
    
    print(paste0("Working on individual ", i, " of ", length(allInd)));
    
    cur_df <- read.table(paste0("derivedAlleleFreqs_snpEff_byInd/", cur_ind_file), sep = "\t", header = T, stringsAsFactors = F);
    
    if (nrow(cur_df) != nrow(sitesInOrder)) {
      # skip this individual if the number of sites is not ok;
      print(paste0("ERROR: number of sites does not match expected number for individual ", cur_ind));
      break;
    }
    
    # subset to only the filtered sites of interest;
    # for this, we will need to use the sites defined in glbl;
    # the paste0() makes the sorting with %in% much faster than sorting by two separate columns in each data.frame without pasting;
    cur_df <- data.frame(gt = cur_df[goodSitesIndex, ]);
    
    if (nrow(cur_df) != nrow(glbl)) {
      # skip this population if the number of individuals is not ok;
      print(paste0("ERROR: not all necessary filtered SNPs are present for population ", cur_pop));
      break;
    }
 
    # combine with the info we need from glbl;
    cur_df <- cbind(glbl, cur_df);
       
    #
    
    # we are not interested in any sites that are missing data;
    cur_df_nonMiss <- cur_df[!(is.na(cur_df$gt)), ];
    
    # split into different categories of interest;
    cur_df_MODIFIER <- cur_df_nonMiss[cur_df_nonMiss$putative_impact == "MODIFIER", ];
    cur_df_LOW <- cur_df_nonMiss[cur_df_nonMiss$putative_impact == "LOW", ];
    cur_df_MODERATE <- cur_df_nonMiss[cur_df_nonMiss$putative_impact == "MODERATE", ];
    cur_df_HIGH <- cur_df_nonMiss[cur_df_nonMiss$putative_impact == "HIGH", ];

    # calculate traits of interest;
    
    out_df <- data.frame(ind = cur_ind);
    
    out_df$MODIFIER_nSitesAll = nrow(cur_df_MODIFIER);
    out_df$MODIFIER_nHomAnc = sum(cur_df_MODIFIER$gt == 0);
    out_df$MODIFIER_nHet = sum(cur_df_MODIFIER$gt == 1);
    out_df$MODIFIER_nHomDer = sum(cur_df_MODIFIER$gt == 2);
    out_df$MODIFIER_propHomAnc = out_df$MODIFIER_nHomAnc / out_df$MODIFIER_nSitesAll;
    out_df$MODIFIER_propHet = out_df$MODIFIER_nHet / out_df$MODIFIER_nSitesAll;
    out_df$MODIFIER_propHomDer = out_df$MODIFIER_nHomDer / out_df$MODIFIER_nSitesAll;

    out_df$LOW_nSitesAll = nrow(cur_df_LOW);
    out_df$LOW_nHomAnc = sum(cur_df_LOW$gt == 0);
    out_df$LOW_nHet = sum(cur_df_LOW$gt == 1);
    out_df$LOW_nHomDer = sum(cur_df_LOW$gt == 2);
    out_df$LOW_propHomAnc = out_df$LOW_nHomAnc / out_df$LOW_nSitesAll;
    out_df$LOW_propHet = out_df$LOW_nHet / out_df$LOW_nSitesAll;
    out_df$LOW_propHomDer = out_df$LOW_nHomDer / out_df$LOW_nSitesAll;

    out_df$MODERATE_nSitesAll = nrow(cur_df_MODERATE);
    out_df$MODERATE_nHomAnc = sum(cur_df_MODERATE$gt == 0);
    out_df$MODERATE_nHet = sum(cur_df_MODERATE$gt == 1);
    out_df$MODERATE_nHomDer = sum(cur_df_MODERATE$gt == 2);
    out_df$MODERATE_propHomAnc = out_df$MODERATE_nHomAnc / out_df$MODERATE_nSitesAll;
    out_df$MODERATE_propHet = out_df$MODERATE_nHet / out_df$MODERATE_nSitesAll;
    out_df$MODERATE_propHomDer = out_df$MODERATE_nHomDer / out_df$MODERATE_nSitesAll;

    out_df$HIGH_nSitesAll = nrow(cur_df_HIGH);
    out_df$HIGH_nHomAnc = sum(cur_df_HIGH$gt == 0);
    out_df$HIGH_nHet = sum(cur_df_HIGH$gt == 1);
    out_df$HIGH_nHomDer = sum(cur_df_HIGH$gt == 2);
    out_df$HIGH_propHomAnc = out_df$HIGH_nHomAnc / out_df$HIGH_nSitesAll;
    out_df$HIGH_propHet = out_df$HIGH_nHet / out_df$HIGH_nSitesAll;
    out_df$HIGH_propHomDer = out_df$HIGH_nHomDer / out_df$HIGH_nSitesAll;
    
    # write output;
    
    write.table(out_df, paste0("outStats_snpEff_byInd/outStats_snpEff_", cur_ind, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }
  
}

#

# combine outfiles into a single file;
# do not print the first line, as it is a header in all files - print the header first as its own separate line;
# then, in the second command, append using ">>";

system(paste0("cat outStats_snpEff_byInd/outStats_snpEff_", allInd[1], ".txt | head -1 > outStats_snpEff_byInd/outStats_allInd_snpEff.txt"));
system(paste0("tail -n +2 -q outStats_snpEff_byInd/outStats_snpEff_*.txt >> outStats_snpEff_byInd/outStats_allInd_snpEff.txt"));
