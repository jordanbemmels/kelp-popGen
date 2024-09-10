##### Script is to generate a new mask file that excludes windows that extremely frequently contain ROHets. These windows are identified by counting the number of times a window contains a ROHet across all individuals, then excluding windows above the 99.99th percentile expected from a Poisson distribution. The new mask can then be used to re-run bcftools roh and re-identify the ROHs, prior to running roh-selfing. #####

# Author use only: GitHub version derived from roh-selfing KL14 script remask_VCF_KL14_getMask_v2.R, 2024/09/10;

#####

### Prerequisites:

# 1) Identification of ROHets ≤10-kbp in size for each individual, previously generated with identify_short_ROHets.R

	# initialRoundROHs_rohet_10kbp.txt;
	
# 2) A file that identifies any 10-kbp windows that were masked from the genome when calculating the ROHets and ROHs, to avoid using these regions in statistics:
# format is four columns: chromosome, start_pos, end_pos, mask [0 for unmasked windows, 1 for masked if the window completely overlaps with masked regions];
# e.g., JARUPZ010000002.1     1 10000    0

	# maskFromRed_10kbp_regions.txt;

#####

# load ROHets;

df <- read.table("initialRoundROHs_rohet_10kbp.txt", sep = "\t", header = T, stringsAsFactors = F);

sampleColumns <- 4:ncol(df);

#########################
#########################

# add mask regions in the data frame for ROHeterozygosity;

mask_10kbp <- read.table("maskFromRed_10kbp_regions.txt", sep = "\t", header = T, stringsAsFactors = F);

df <- merge(df, mask_10kbp, sort = F);

#########################
#########################

##### Summarize by individual #####

# user may need to edit this section if samples are named with a different convention;
byInd <- data.frame(ind = colnames(df)[sampleColumns]);
byInd$ID <- gsub("KL..\\.", "", byInd$ind);
byInd$ID <- gsub("\\.", "-", byInd$ID);

byInd$ROHet_windows <- apply(df[ , sampleColumns], 2, sum);

# define individuals to exclude;
# user may need to edit this section if different samples are to be excluded;
# here, these samples were excluded because the ROHs and ROHets were determined from bcftools roh command run on <4 individuals using rangewide allele frequencies - inaccurate and should not be included;
toExcludeIDs <- byInd$ID[substr(byInd$ID, 1, 5) %in% c("NL-HY") | substr(byInd$ID, 1, 8) %in% c("NL-CR-03", "NL-CR-04", "NL-CR-05", "NL-CR-06", "NL-CR-07")];

sampleColumns_good <- sampleColumns[!(sampleColumns %in% (3 + which(byInd$ID %in% toExcludeIDs)))];
colnames(df)[sampleColumns_good];

##### Summarize by window #####

byWindow <- df[ , c("chr", "start", "end", "mask")];

byWindow$nInd_ROHet <- apply(df[ , sampleColumns_good], 1, sum);

###

# mask should be the windows that are overrepresented in the number of ROHeterozygosity;
# assume a Poisson distribution;

# what quantile should we use? check how many windows above a certain quantile cutoff that we expect by chance;
nrow(byWindow[byWindow$mask != 1, ]) * 0.001; # expect 38 windows; 
nrow(byWindow[byWindow$mask != 1, ]) * 0.0001; # expect 3.8 windows - seems a reasonable cutoff;

qpois(0.9999, mean(byWindow$nInd_ROHet[byWindow$mask != 1])); #10;

# thus, we will mask any regions in which >10 individuals are in a short ROHet (with short ROHet defined as a ROHet ≤ 10kbp);
# note: why are we using >10 instead of >=10? because according to ?qpois(), "The quantile is right continuous: qpois(p, lambda) is the smallest integer x such that P(X<=x) >= p." So if x=10 then the probability that the integer is 10 or less includes p itself (">=" p), but we want to exclude the probability of p itself (0.9999), so we need to go up to the next highest integer using >10;

# how many windows is this, and what proportion?;
nrow(byWindow[byWindow$mask != 1 & byWindow$nInd_ROHet > 10, ]); # 1576;
nrow(byWindow[byWindow$mask != 1 & byWindow$nInd_ROHet > 10, ]) / nrow(byWindow[byWindow$mask != 1, ]); # 0.0414;

# how many short ROHet will actually be removed total (approximately, the HMM will recalculate ROHs);
sum(byWindow$nInd_ROHet[byWindow$mask != 1 & byWindow$nInd_ROHet > 10]); # 30,084;

# visualize;
hist(byWindow$nInd_ROHet[byWindow$mask != 1], breaks = 100, xlab = "Number of individuals for which a window is in an ROHeterozygosity", ylab = "Number of windows", main = "Red line = 0.9999 quantile in Poisson distribution with same mean");
abline(v = 10, col = "red");

#########################
#########################

##### Create a new mask and save the output #####

byWindow$newMask <- 0;
byWindow$newMask[byWindow$mask != 1 & byWindow$nInd_ROHet > 10] <- 1;

# to avoid writing exponential notation, we need set scipen to a very high value;
options(scipen = 999);

write.table(byWindow[byWindow$newMask == 1, c("chr", "start", "end")], file = "remask_VCF_mask_KL14_v2.bed", quote = F, sep = "\t", col.names = F, row.names = F);
