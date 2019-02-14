# randInte_IntersectionMultiValidRefs.R
#Script gets the intersection of the 2 valid het Ref overlap.rda files provided. 
#
#ref1 - e.g. GWAS -SCZ, PD and
#ref2 - SPIDEX, Annotation - transcript specific


#R --vanilla --slave < ./randInte_IntersectionMultiValidRefs.R <ref1Overlap.rda> <ref2Overlap.rda> <out.rda>
#where
#ref1Overlap.rda = validhets/queryset SNPs overlap with ref1 - e.g. GWAS -SCZ, PD and
#ref2Overlap.rda = validhets/queryset SNPs overlap with ref2 - e.g. SPIDEX, Annotation - transcript specific
#out.rda        = output to be saved with this file name.

#Output saved in <BASEDIR>files/overlapFiles/

#rm(list=ls())
ref1.rda <- as.character( commandArgs()[4] )
ref2.rda <- as.character( commandArgs()[5] )
out.rda <- as.character( commandArgs()[6] ) 
#refData  <- as.character( commandArgs()[7] )

BASEDIR <- "/home/kdsa/bootstrap/"
source(paste0(BASEDIR, "scripts/randInteMultiRefs_functions.R"))
source(paste0(BASEDIR, "scripts/process_ASEdatasets_functions.R"))

outdir <- paste0(BASEDIR, "files/overlapFiles/")
setwd(outdir)

#PD, spidex
#ref1.rda <- "PUTM_validhets_PD2018Meta5Ex23andMe_Overlap.rda"
#ref2.rda <- "PUTM_validhets_SPIDEX_Overlap.rda"
#out.rda <- "PUTM_validhets_PD2018Meta5Ex23andMe_SPIDEX_Overlap.rda"

cat("\nRef1 = ", ref1.rda)
cat("\nRef2 = ", ref2.rda)
cat("\nOutput file = ", out.rda)
cat("\n in ", getwd(), "\n")
#calling function to get the intersection - for diff ref1s
getIntersectionValRefs(ref1.rda, ref2.rda, out.rda )

#to check
#load("PUTM_validhets_PD2018Meta5Ex23andMe_SPIDEX_Overlap.rda")
#head(valid.ref)

#checking the distribution of Avg reads in the ASEs vs non ASEs
bwidthlist <- 2
generatePlot_AvgReadDistribution4Bin(out.rda ,bwidthlist, "F" ) # single ref <- "F"
