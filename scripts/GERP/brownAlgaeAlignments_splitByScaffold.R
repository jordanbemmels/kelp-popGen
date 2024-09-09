##### The script takes multiple, separate fasta-format alignment files (aligned to the focal species' reference genome) for each brown algae genome and splits them into each scaffold. It then combines the alignments for each scaffold across species, with the end result being one alignment file per scaffold that contains all brown algae species. #####

# Author use only: GitHub version derived from KL14 script 03_splitByScaffold.R, 2024/09/09;

# Code is rewritten with heavy modifications from a script by Becky S Taylor:
# https://github.com/BeckySTaylor/Phylogenomic_Analyses/blob/main/GERP_running

# User sets these manually each time;

pathToRetainedSpeciesList <- "/path/to/retained/species/list.txt";
pathToAutosomalScaffoldsList <- "/path/to/autosomal/scaffolds/list.txt";

#####

system("mkdir byScaffold", intern = F);
system("mkdir byScaffold/GERP_formatted", intern = F);

# load a list of all brown algage species that have been individually aligned to the focal species' reference genome;
# the retained species list is a text file, with one species per line;
# this assumes that the alignments for each brown algae species are in the current working directory and named like "species01_filtered_sorted.fasta", "species02_filtered_sorted.fasta", etc.;
allSpp <- read.table(pathToRetainedSpeciesList, sep = "\t", header = F, stringsAsFactors = F);

# load a list of the autosomal scaffolds from the focal species' genome, one scaffold per line;
autoScf <- read.table(pathToAutosomalScaffoldsList, header = F, sep = "\t", stringsAsFactors = F);

# for each species, split the alignments into separate scaffolds;

for (i in 1:nrow(allSpp)) {
  
  cur_sp <- allSpp$V1[i];
    
  # split into separate files, one for each scaffold;
  cmd_csplit <- paste0("csplit -s -z ", cur_sp, "_filtered_sorted.fasta '/>/' '{*}'");
  system(cmd_csplit, intern = F);
  
  # check the first line of each scaffold file to determine which scaffold it is in (this cannot be predicted 1:1 based on the file name);
  all_xx <- list.files("./", "xx");
  for (j in 1:length(all_xx)) {
    con <- file(all_xx[j]);
    first_line <- readLines(con, n=1);
    close(con);
    #
    cur_scf <- gsub(">", "", first_line);
    #
    # change first line to contain the species name - own modification of code;
    # also write to a new outfile named after the species and scaffold;
    cmd_sed <- paste0("sed '1 s/^.*$/>", cur_sp, "/' ", all_xx[j], " > byScaffold/", cur_sp, "_filtered_sorted_", cur_scf, ".fasta");
    system(cmd_sed, intern = F);
    #
    # remove the corresponding xx file;
    cmd_rm_xx <- paste0("rm ", all_xx[j]);
    system(cmd_rm_xx, intern = F);
  }
  
}

# once each scaffold file has been renamed and processed for all individuals, make one fasta file for each scaffold - only process the autosomal scaffolds of interest;
for (i in 1:nrow(autoScf)) {
  cmd_cat <- paste0("cat byScaffold/*_", autoScf$V1[i], ".fasta > byScaffold/GERP_formatted/", autoScf$V1[i], "_GERP.fasta");
  system(cmd_cat, intern = F);
  #
  # reformat fasta file for this program as each individual sequence needs to be on one line, not multiple lines;
  # note that r"()" allows both single and double quotes to be used in a string, which is necessary here for awk;
  # I modified the command because nohup.out wasn't actually automatically generated when running through R, so instead print to an outfile that contains the actual name we want, avoiding the nohup and move commands;
  cmd_nohup <- paste0(r"(awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' byScaffold/GERP_formatted/)", autoScf$V1[i], "_GERP.fasta > byScaffold/GERP_formatted/", autoScf$V1[i], "_GERP_formatted.mfa");
  system(cmd_nohup, intern = F);
  #
  # remove the first line that is a blank;
  cmd_sed_rm <- paste0("sed '1d' -i byScaffold/GERP_formatted/", autoScf$V1[i], "_GERP_formatted.mfa");
  system(cmd_sed_rm, intern = F);
  #
  # remove the file that was a fasta with sequences on multiple lines;
  cmd_rm_fasta <- paste0("rm byScaffold/GERP_formatted/", autoScf$V1[i], "_GERP.fasta");
  system(cmd_rm_fasta, intern = F);
}

# remove all individual scaffold files (that existed prior to merging into a single file for all individuals);
cmd_rm_indScf <- paste0("rm byScaffold/*.fasta");
system(cmd_rm_indScf, intern = F);
