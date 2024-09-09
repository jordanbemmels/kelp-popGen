##### Script is to identify short Runs of Heterozygosity (ROHets), not to be confused with Runs of Homozygosity (ROHs). The ROHets are short regions of heterozygosity surrounded by ROHs. The script will calculate ROHets that are up to 1 kbp, 10 kbp, and 100 kbp in length. #####

# Author use only: GitHub version derived from roh-selfing KL14 script initialRoundROHs_KL14_summarize_v2.R, 2024/09/09;

#####

### Prerequisites:

# 1) User has already identified ROHs with bcftools for each population, with output stored in the initialRoundROHs directory. This can be done as follows, for example for population pop_01 from a vcf file for that population:

	# bcftools roh -O r -o initialRoundROHs/pop_01_rohoutput.txt pop_01.vcf

###

# User sets these manually each time;

pathToAutosomalScaffoldsList <- "/path/to/autosomal/scaffolds/list.txt"; # note that it must have the length of the scaffolds in the third column;

#####

rohFiles <- list.files("initialRoundROHs");

scaffolds <- read.table(pathToAutosomalScaffoldsList, sep = "\t", header = F, stringsAsFactors = F);

# set up a dataframe in 10-kbp blocks;

pos_start <- seq(1, scaffolds$V3[1], 10000);
pos_end <- pos_start + 9999;
chr <- rep(scaffolds$V1[1], length(pos_start));

for (i in 2:nrow(scaffolds)) {
  
  cur_pos_start <- seq(1, scaffolds$V3[i], 10000);
  cur_pos_end <- cur_pos_start + 9999;
  cur_chr <- rep(scaffolds$V1[i], length(cur_pos_start));
  
  pos_start <- c(pos_start, cur_pos_start);
  pos_end <- c(pos_end, cur_pos_end);
  chr <- c(chr, cur_chr);
  
}

df <- data.frame(chr=chr, start=pos_start, end=pos_end);
# divide into different data frames to store results of different run of HETEROzygosity lengths;
df_1kbp <- df;
df_10kbp <- df;
df_100kbp <- df;

# process each initial round of ROHs file;

for (i in 1:length(rohFiles)) {
  
  print(paste0("Working on file ", i , " of ", length(rohFiles)));
  
  roh <- read.table(paste0("initialRoundROHs/", rohFiles[i]), sep = "\t", header = F, stringsAsFactors = F);
  colnames(roh) <- c("RG", "Sample", "Chromosome", "Start", "End", "Length", "Number_of_markers", "Quality");
  
  cur_samples <- unique(roh$Sample);
  
  for (j in 1:length(cur_samples)) {

    print(paste0("Working on sample ", j , " of ", length(cur_samples)));
        
    cur_roh <- roh[roh$Sample == cur_samples[j], ];
    df_1kbp[ , ncol(df_1kbp) + 1] <- 0;
    df_10kbp[ , ncol(df_10kbp) + 1] <- 0;
    df_100kbp[ , ncol(df_100kbp) + 1] <- 0;
    colnames(df_1kbp)[ncol(df_1kbp)] <- cur_samples[j];
    colnames(df_10kbp)[ncol(df_10kbp)] <- cur_samples[j];
    colnames(df_100kbp)[ncol(df_100kbp)] <- cur_samples[j];
    
    for (k in 1:nrow(scaffolds)) {
      
      cur_scf <- scaffolds$V1[k];
      cur_roh_scf <- cur_roh[cur_roh$Chromosome == cur_scf, ];
      
      if (nrow(cur_roh_scf) > 0) {
        
        cur_rohet <- data.frame(Chromosome = cur_scf, Start = 1, End = cur_roh_scf$Start[1] - 1);
        
        if (nrow(cur_roh_scf) > 2) {
          
          for (m in 2:nrow(cur_roh_scf)) {
            cur_rohet <- rbind(cur_rohet, c(cur_scf, cur_roh_scf$End[m - 1] + 1, cur_roh_scf$Start[m] - 1));
          }
          
        }
        
        cur_rohet <- rbind(cur_rohet, c(cur_scf, min(cur_roh_scf$End[nrow(cur_roh_scf)] + 1, scaffolds$V3[k]), scaffolds$V3[k]));
 
        cur_rohet$Start <- as.numeric(cur_rohet$Start);
        cur_rohet$End <- as.numeric(cur_rohet$End);
               
        cur_rohet$Length <- cur_rohet$End - cur_rohet$Start + 1;
        cur_rohet <- cur_rohet[cur_rohet$Length > 0, ];
        
        # cur_rohet is ready for processing;
        
        for (n in 1:nrow(cur_rohet)) {
          
          # print the run of HETERozygosity that overlaps even slightly with the interval;
          
          if (cur_rohet$Length[n] <= 1000) {
            df_1kbp[df_1kbp$chr == cur_rohet$Chromosome[n] & !(df_1kbp$start < cur_rohet$Start[n] - 9999) & !(df_1kbp$end > cur_rohet$End[n] + 9999), ncol(df_1kbp)] <- 1;
          }
          
          if (cur_rohet$Length[n] <= 10000) {
            df_10kbp[df_10kbp$chr == cur_rohet$Chromosome[n] & !(df_10kbp$start < cur_rohet$Start[n] - 9999) & !(df_10kbp$end > cur_rohet$End[n] + 9999), ncol(df_10kbp)] <- 1;
          }
          
          if (cur_rohet$Length[n] <= 100000) {
            df_100kbp[df_100kbp$chr == cur_rohet$Chromosome[n] & !(df_100kbp$start < cur_rohet$Start[n] - 9999) & !(df_100kbp$end > cur_rohet$End[n] + 9999), ncol(df_100kbp)] <- 1;
          }

        }
        
      }
      
    }  

  }
  
}

write.table(df_1kbp, "initialRoundROHs_rohet_1kbp.txt", sep = "\t", row.names = F, col.names = T, quote = F);

write.table(df_10kbp, "initialRoundROHs_rohet_10kbp.txt", sep = "\t", row.names = F, col.names = T, quote = F);

write.table(df_100kbp, "initialRoundROHs_rohet_100kbp.txt", sep = "\t", row.names = F, col.names = T, quote = F);
