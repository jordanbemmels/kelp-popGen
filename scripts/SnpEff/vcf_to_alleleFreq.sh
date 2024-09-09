#!/bin/bash
#Code is to parse a vcf file with annotations added by SnpEff,\
#and print the allele frequencies (across a subset of samples\
#defined by the sample list, provided as an input file) plus\
#some basic info about the annotation

#Usage: bash this_script.sh vcfFile sampleFile outFile

vcfFile=$1;
sampleFile=$2;
outFile=$3;

#####

bcftools view -H $1 \
	-S $2 \
	| awk 'BEGIN {print "chr", "pos", "homRef", "het", "homAlt", "missing", "total", "freqRef", "freqAlt", "annotation", "putative_impact", "gene_name"} \
	{homRef=0; het=0; homAlt=0; missing=0; z=match($8 , /ANN=[^;]*/); \
	for(i=10; i<=NF; i++) {homRef += gsub(/0\/0/,1,$i); \
	het += gsub(/0\/1/,1,$i); \
	homAlt += gsub(/1\/1/,1,$i); \
	missing += gsub(/\.\/\./,1,$i)}; \
	if(z) { ann=substr($8, RSTART, RLENGTH)}; \
	split(ann,a,"[|]"); \
	print $1, $2, homRef, het, homAlt, missing, homRef+het+homAlt+missing, (2*homRef+het)/(2*(homRef+het+homAlt)), (het+2*homAlt)/(2*(homRef+het+homAlt)), a[2], a[3], a[4]}' \
	> $3
