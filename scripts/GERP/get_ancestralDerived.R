##### This script determines which allele is ancestral and derived. The script is run after running the modified gerpcol script from Taylor et al. 2024, as described in the manuscript. #####

# Author use only: GitHub version derived from KL14 script get_ancestralDerived_GERP_KL14.R, 2024/09/09;

# This example is to determine the ancestral/derived allele for bull kelp (Nereocystis). If the focal species in instead giant kelp (Macrocystis), the user needs to change "Macrpy" to "Nerelu" throughout this script, so that the outgroup will be set to Nereocystis luetkeana.;

# User sets these manually each time;

pathToAutosomalScaffoldsList <- "/path/to/autosomal/scaffolds/list.txt";
pathToRefAltDefinitions <- "path/to/ref/alt/definitions.txt";
outputFileName <- "outfile.txt";

#####

# load a list of the autosomal scaffolds from the focal species' genome, one scaffold per line;
autoScf <- read.table(pathToAutosomalScaffoldsList, header = F, sep = "\t", stringsAsFactors = F);

# read in all the variable sites in the genome (these correspond to a VCF file for which you are interested in defining the ancestral and derived alleles), with the reference and alternative alleles printed;
# the format for the reference and alternative allele definitions is a text file with four columns [chromosome, position, ref-allele, alt-allele], e.g., the first line may be:
	# JARUPZ010000002.1	16303	G	A

df <- read.table(pathToRefAltDefinitions, header = F, sep = "\t", stringsAsFactors = F);
colnames(df) <- c("chr", "pos", "ref", "alt");

#

# combine gerpcol output for all scaffolds with the ref and alt allele definitions for each genomic site of interest;

for (i in 1:nrow(autoScf)) {
  
  # counter;
  print(paste0("Reading in gerpcol for scaffold ", autoScf$V1[i]));
  
  # read in the gerpcol output - relevant lines only for each scaffold;
  # the gerpcol output has been modified to contain only the lines corresponding to the genomic sites of interest in the VCF file for which ancestral and derived alleles are required, and corresponding to the same lines as in pathToRefAltDefinitions above;
  # there is one gerpcol output file per scaffold, in the directory gerpcol_relevantLines_byScaffold;
  gerpcol <- read.table(paste0("gerpcol_relevantLines_byScaffold/gerpcol_", autoScf$V1[i], ".txt"), header = F, sep = "\t", stringsAsFactors = F);
  colnames(gerpcol) <- c("pos", "N", "S", "Macrpy", "Saccja", "Lamidi");
  
  # change the gerpcol pos to be 1-indexed instead of 0-indexed, to match the format of the original VCF indexing and df;
  gerpcol$pos <- gerpcol$pos + 1;
  
  # add chromosome to gerpcol;
  gerpcol$chr <- autoScf$V1[i];
 
  # subset df to only the scaffold of interest;
  cur_df <- df[df$chr == autoScf$V1[i], ];
  
  # merge with gerpcol;
  # UPDATE 2024/04/23 added "sort = F" to ensure that the order of alleles is not changed; 
  cur_df <- merge(cur_df, gerpcol, all = T, sort = F);
  
  # combine into output;
  if (i == 1) {
    out <- cur_df;
  } else {
    out <- rbind(out, cur_df);
  }
   
}

#####

# define the ancestral and derived allele - in order for an allele to be derived, it must not be present in any of the three outgroups (but missing data is allowed, i.e., we do not need data for all three outgroups);
# if neither the ref nor the alt is non-present in any of the three outgroups, then we cannot define ancestral and derived alleles;

out$N <- as.numeric(out$N);
out$S <- as.numeric(out$S);

out$n_outgroup_alleles <- apply(out[ , c("Macrpy", "Saccja", "Lamidi")], 1, function(x) {length(unique(x[!is.na(x) & x != "N"]))});

out$unique_outgroup_allele <- NA;
out$unique_outgroup_allele[out$n_outgroup_alleles == 1] <- apply(out[out$n_outgroup_alleles == 1, c("Macrpy", "Saccja", "Lamidi")], 1, function(x) {unique(x[!is.na(x) & x != "N"])});

out$ancestral <- NA;
out$derived <- NA;

out$ancestral[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)] <- out$ref[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)];
out$ancestral[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)] <- out$alt[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)];

out$derived[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)] <- out$alt[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)];
out$derived[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)] <- out$ref[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)];

#####

write.table(out, outputFileName, sep = "\t", row.names = F, col.names = T, quote = F);
