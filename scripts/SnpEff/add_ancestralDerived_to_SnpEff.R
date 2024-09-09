##### This script adds the ancestral vs. derived definitions (previously determined during GERP analyses) to the partially-processed SnpEff annotations. For processing efficiency, the input should be split by scaffold (see prerequisites section), but the output is combined for all scaffolds. #####

# Author use only: GitHub version derived from KL14 script add_ancestralDerived_to_SnpEff_KL14_KRI.R, 2024/09/09;

### Prerequisites:

# 1) Custom ref-alt allele frequency file with SnpEff annotations (alleleFreq.txt) need to split by chromosome and each split file placed within the byScaffold/ directory. A file for the whole genome was generated previously with vcf_to_alleleFreq.sh. It can be split by scaffold, for example if the first scaffold is chr_1:

	# awk '$1 == \"chr_1\" {print $0}' alleleFreq.txt > byScaffold/chr_1_alleleFreq.txt

###

# User sets these manually each time;

pathToAutosomalScaffoldsList <- "/path/to/autosomal/scaffolds/list.txt";
pathToAncestralDerivedDefinitions <- "/path/to/ancestralDerived.txt"; # previously generated during GERP analysis;

#####

system("mkdir temp", intern = F);

autoScf <- read.table(pathToAutosomalScaffoldsList, header = F, sep = "\t", stringsAsFactors = F);

# load ancestral-derived file;

anc <- read.table(pathToAncestralDerivedDefinitions, sep = "\t", header = T, stringsAsFactors = F);

# process each split-scaffold SnpEff file;

for (i in 1:nrow(autoScf)) {
  
  print(paste0("Working on scaffold ", i, " of ", nrow(autoScf)));
  
  snpEff <- read.table(paste0("byScaffold/", autoScf$V1[i], "_alleleFreq.txt"), sep = " ", header = F, stringsAsFactors = F);
  colnames(snpEff) <- c("chr", "pos", "homRef", "het", "homAlt", "missing", "total", "freqRef", "freqAlt", "annotation", "putative_impact", "gene_name");
  
  # reduce anc to the current scaffold;
  anc_cur_scf <- anc[anc$chr == autoScf$V1[i], ];
  
  # we don't need to keep info about GERP score anymore, because this analysis does not use GERP and only wants to ancestral/derived alleles;
  anc_cur_scf <- anc_cur_scf[ , c("chr", "pos", "ref", "alt", "ancestral", "derived")];
  
  # merge the snpEff and ancestral/derived sets;
  # use all = F because we are only interested in site that BOTH have the ancestral allele defined, and have a SnpEff annotation;
  # there will be many sites without an ancestral allele defined, because there were not enough alignments from the outgroup species when defining ancestral alleles;
  # we there expect the number of sites to be reduced substantially;
  cur_merged_df <- merge(anc_cur_scf, snpEff, all = F);
  
  # sort by the position;
  cur_merged_df <- cur_merged_df[order(cur_merged_df$pos), ];
  
  # combine with an output dataframe;
  if (i == 1) {
    df_out <- cur_merged_df;
  } else {
    df_out <- rbind(df_out, cur_merged_df);
  }
  
}

# in some cases, the ancestral/derived is actually missing - this occurs when the unique outgroup allele was *different* from either of the two alleles for the SNP that are segregating in the focal species;

df_out <- df_out[!(is.na(df_out$ancestral)), ];


# add a couple other columns of useful info to the output dataframe;

df_out$ref_is_anc <- 0;
df_out$ref_is_anc[df_out$ref == df_out$ancestral] <- 1;

unique(df_out[df_out$ref_is_anc == 1, c("ref", "ancestral")]);
unique(df_out[df_out$ref_is_anc == 0, c("ref", "ancestral")]);

df_out$freqDer <- NA;
df_out$freqDer[df_out$ref_is_anc == 1] <- df_out$freqAlt[df_out$ref_is_anc == 1];
df_out$freqDer[df_out$ref_is_anc == 0] <- df_out$freqRef[df_out$ref_is_anc == 0];

# subset to the columns we need to retain;

df_out <- df_out[ , c("chr", "pos", "ref", "alt", "ancestral", "derived", "freqAlt", "freqDer", "annotation", "putative_impact", "gene_name")];

# write output;

write.table(df_out, "alleleFreq_ancestralDerived.txt", sep = "\t", row.names = F, col.names = T, quote = F);

