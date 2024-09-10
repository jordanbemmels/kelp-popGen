##### Script is to calculate realized genetic load in each individual. This script processes the GERP analysis, and divides the sites into evolutionarily labile vs. evolutionarily conserved sites. #####

# Author use only: GitHub version derived from Part II only (individual-level analysis) from version 5 of KL14 GERP script analyzeOutput_GERP_KL14_ancestralGERP_v5.R, 2024/09/10;

#####

### Prerequisites:

# 1) Derived allele frequencies separately for each population, and separately for each individual, all calculated previously using the GERP script get_derivedAlleleFreqs.R:

	# derivedAlleleFreqs_byPop/derivedAlleleFreqs_[...].txt; # results for each population [...];
	# derivedAlleleFreqs_byInd/derivedAlleleFreqs_[...].txt; # results for each individual [...];

# 2) List of filtered sites to use for calculating realized load, previously calculated using the purging_drift_GERP.R script:

	# global_filteredSitesList_GERP.txt;
	# global_filteredSitesList_min2_GERP.txt;
	
#########################
#########################

library(foreach);
library(doParallel);

##############################
##############################

##### PROCESS BY INDIVIDUAL #####

library(foreach);
library(doParallel);

# re-load our lists of good sites;
glbl <- read.table("global_filteredSitesList_GERP.txt", header = T, sep = "\t", stringsAsFactors = F);
glbl_min2 <- read.table("global_filteredSitesList_min2_GERP.txt", header = T, sep = "\t", stringsAsFactors = F);

# re-load our list of all sites, so that we know which line is which site in the individual files - this was not outputted, in order to save space;
# just load the first population as an example, here NL-BA-01 - we don't need all populations;
sitesInOrder <- read.table("derivedAlleleFreqs_byPop/derivedAlleleFreqs_NL-BA-01.txt", sep = "\t", header = T, stringsAsFactors = F);
sitesInOrder <- sitesInOrder[ , c("chr", "pos")];

goodSitesIndex_min0 <- which(paste0(sitesInOrder$chr, "_", sitesInOrder$pos) %in% paste0(glbl$chr, "_", glbl$pos));
length(goodSitesIndex_min0) == nrow(glbl);

goodSitesIndex_min2 <- which(paste0(sitesInOrder$chr, "_", sitesInOrder$pos) %in% paste0(glbl_min2$chr, "_", glbl_min2$pos));
length(goodSitesIndex_min2) == nrow(glbl_min2);

# define all individuals;

allIndFiles <- list.files("derivedAlleleFreqs_byInd");
allInd <- gsub("derivedAlleleFreqs_", "", allIndFiles);
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
    
    print(paste0("Working on population ", i, " of ", length(allInd)));
    
    cur_df <- read.table(paste0("derivedAlleleFreqs_byInd/", cur_ind_file), sep = "\t", header = T, stringsAsFactors = F);

    if (nrow(cur_df) != nrow(sitesInOrder)) {
      # skip this individual if the number of sites is not ok;
      print(paste0("ERROR: number of sites does not match expected number for individual ", cur_ind));
      break;
    }
    
    # subset to only the filtered sites of interest;
    # for this, we will need to use the sites defined in glbl;
    # the paste0() makes the sorting with %in% much faster than sorting by two separate columns in each data.frame without pasting;
    cur_df_min0 <- data.frame(gt = cur_df[goodSitesIndex_min0, ]);
    cur_df_min2 <- data.frame(gt = cur_df[goodSitesIndex_min2, ]);
    
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
     
    # combine with the info we need from glbl;
    cur_df_min0 <- cbind(glbl, cur_df_min0);
    cur_df_min2 <- cbind(glbl_min2, cur_df_min2);
           
    #
    
    # we are not interested in any sites that are missing data;
    cur_df_min0_nonMiss <- cur_df_min0[!(is.na(cur_df_min0$gt)), ];
    cur_df_min2_nonMiss <- cur_df_min2[!(is.na(cur_df_min2$gt)), ];
    
    # split into different categories of interest;
    cur_df_min0_LABILE <- cur_df_min0_nonMiss[cur_df_min0_nonMiss$N > 0.5 & cur_df_min0_nonMiss$S <= 0.5, ];
    cur_df_min0_CONSERVED <- cur_df_min0_nonMiss[cur_df_min0_nonMiss$N > 0.5 & cur_df_min0_nonMiss$S > 0.5, ];

    cur_df_min2_LABILE <- cur_df_min2_nonMiss[cur_df_min2_nonMiss$N > 0.5 & cur_df_min2_nonMiss$S <= 0.5, ];
    cur_df_min2_CONSERVED <- cur_df_min2_nonMiss[cur_df_min2_nonMiss$N > 0.5 & cur_df_min2_nonMiss$S > 0.5, ];
        
    # calculate traits of interest;
    
    out_df_min0 <- data.frame(ind = cur_ind);
        
    out_df_min0$LABILE_nSitesAll = nrow(cur_df_min0_LABILE);
    out_df_min0$LABILE_nHomAnc = sum(cur_df_min0_LABILE$gt == 0);
    out_df_min0$LABILE_nHet = sum(cur_df_min0_LABILE$gt == 1);
    out_df_min0$LABILE_nHomDer = sum(cur_df_min0_LABILE$gt == 2);
    out_df_min0$LABILE_propHomAnc = out_df_min0$LABILE_nHomAnc / out_df_min0$LABILE_nSitesAll;
    out_df_min0$LABILE_propHet = out_df_min0$LABILE_nHet / out_df_min0$LABILE_nSitesAll;
    out_df_min0$LABILE_propHomDer = out_df_min0$LABILE_nHomDer / out_df_min0$LABILE_nSitesAll;

    out_df_min0$CONSERVED_nSitesAll = nrow(cur_df_min0_CONSERVED);
    out_df_min0$CONSERVED_nHomAnc = sum(cur_df_min0_CONSERVED$gt == 0);
    out_df_min0$CONSERVED_nHet = sum(cur_df_min0_CONSERVED$gt == 1);
    out_df_min0$CONSERVED_nHomDer = sum(cur_df_min0_CONSERVED$gt == 2);
    out_df_min0$CONSERVED_propHomAnc = out_df_min0$CONSERVED_nHomAnc / out_df_min0$CONSERVED_nSitesAll;
    out_df_min0$CONSERVED_propHet = out_df_min0$CONSERVED_nHet / out_df_min0$CONSERVED_nSitesAll;
    out_df_min0$CONSERVED_propHomDer = out_df_min0$CONSERVED_nHomDer / out_df_min0$CONSERVED_nSitesAll;

    #
    
    out_df_min2 <- data.frame(ind = cur_ind);
    
    out_df_min2$LABILE_nSitesAll = nrow(cur_df_min2_LABILE);
    out_df_min2$LABILE_nHomAnc = sum(cur_df_min2_LABILE$gt == 0);
    out_df_min2$LABILE_nHet = sum(cur_df_min2_LABILE$gt == 1);
    out_df_min2$LABILE_nHomDer = sum(cur_df_min2_LABILE$gt == 2);
    out_df_min2$LABILE_propHomAnc = out_df_min2$LABILE_nHomAnc / out_df_min2$LABILE_nSitesAll;
    out_df_min2$LABILE_propHet = out_df_min2$LABILE_nHet / out_df_min2$LABILE_nSitesAll;
    out_df_min2$LABILE_propHomDer = out_df_min2$LABILE_nHomDer / out_df_min2$LABILE_nSitesAll;
    
    out_df_min2$CONSERVED_nSitesAll = nrow(cur_df_min2_CONSERVED);
    out_df_min2$CONSERVED_nHomAnc = sum(cur_df_min2_CONSERVED$gt == 0);
    out_df_min2$CONSERVED_nHet = sum(cur_df_min2_CONSERVED$gt == 1);
    out_df_min2$CONSERVED_nHomDer = sum(cur_df_min2_CONSERVED$gt == 2);
    out_df_min2$CONSERVED_propHomAnc = out_df_min2$CONSERVED_nHomAnc / out_df_min2$CONSERVED_nSitesAll;
    out_df_min2$CONSERVED_propHet = out_df_min2$CONSERVED_nHet / out_df_min2$CONSERVED_nSitesAll;
    out_df_min2$CONSERVED_propHomDer = out_df_min2$CONSERVED_nHomDer / out_df_min2$CONSERVED_nSitesAll;
    
    # write output;
    
    write.table(out_df_min0, paste0("outStats_GERP_byInd/outStats_min0_", cur_ind, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    write.table(out_df_min2, paste0("outStats_GERP_byInd/outStats_min2_", cur_ind, ".txt"), sep = "\t", row.names = F, col.names = T, quote = F);
    
  }
  
}

#

# combine outfiles into a single file;
# do not print the first line, as it is a header in all files - print the header first as its own separate line;
# then, in the second command, append using ">>";

system(paste0("cat outStats_GERP_byInd/outStats_min0_", allInd[1], ".txt | head -1 > outStats_GERP_byInd/outStats_allInd_min0.txt"));
system(paste0("tail -n +2 -q outStats_GERP_byInd/outStats_min0_*.txt >> outStats_GERP_byInd/outStats_allInd_min0.txt"));

system(paste0("cat outStats_GERP_byInd/outStats_min2_", allInd[1], ".txt | head -1 > outStats_GERP_byInd/outStats_allInd_min2.txt"));
system(paste0("tail -n +2 -q outStats_GERP_byInd/outStats_min2_*.txt >> outStats_GERP_byInd/outStats_allInd_min2.txt"));
