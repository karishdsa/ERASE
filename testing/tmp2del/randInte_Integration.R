#script : randInte_Integration.R
#Run integrate() to integrate the transformed zscores of ref1 and ref2 for a set of alphas
#                 alpha = measure of confidence in the ref1 dataset ( between 0-1)
##                 no. of iterations for the random pairwise integration step =10000
#
#R --vanilla --slave < randInte_Integration.R tiss ref1 ref2
#                  tiss - the tissue ( e.g PUTM )
#                  ref1 - the GWAS dataset e.g SCZ
#                  ref2 - the 2nd ref dataset e.g SPIDEX, Annotation

#  Output generated in : <BASEDIR>results/<tiss>_<ref1>_<ref2>
#  with prefix : <tiss>_ASEsNonASEs_<ref1><ref2>_Randomisation_nIter<nIterations>_alpha<alpha>


tiss <- as.character( commandArgs()[4] )
ref1 <- as.character( commandArgs()[5] )
ref2 <- as.character( commandArgs()[6] )
#rm(list=ls())

# ref1 <- "SCZ2018";  ref2 <- "SPIDEX"; tiss <- "PUTM"

BASEDIR <- "/home/kdsa/bootstrap/"
RANDPAIRINT.REPEAT <-TRUE  # random pairwise integration = allowing repeat. Should not make a difference as selecting 1 element at a time

source(paste0(BASEDIR, "scripts/randInteMultiRefs_functions.R"))
cat("\nStart ", date(), "\n")

inp_outdir <- paste0(BASEDIR, "results/", tiss, "_", ref1, "_", ref2)
setwd( inp_outdir)
refs <- paste0(ref1, "_", ref2)
randFilesPattern <- "*obsRandomValsMean_zscoreNmlDistr.rda"
randFiles <- list.files(pattern = randFilesPattern)

alphaVals <- seq(0,1, by=0.1)
nIterations <- 10000


cat("\nTissue :", tiss, "\n")
cat("\nRef1 :", ref1, "\n")
cat("\nRef2 :", ref2, "\n")

cat("\nNo. of iterations :", nIterations)
cat("\nAlpha Values :", alphaVals, "\n")

#check if 2 files found are for the 2 refs       
if ( (length(randFiles) != 2) | ( length(which( grepl (refs, randFiles))) != 2))  {
  cat("\n Randomisation files expected for 2 references.", length(randFiles)," files found\n")
} else {
  #run integration for all alphas
  ret <- lapply(alphaVals, function(alpha){
          #alpha <- 0.5
          cat("\n Integration for alpha= ", alpha, "start..\n")
          val <-runIntegration4alpha ( randFiles[1],randFiles[2], alpha, nIterations, ref1, ref2, tiss) 
          #returns the output rda file name
          cat("\n Integration for alpha= ", alpha, "complete")
          return(val)
        })

}

cat("\nEnd ", date(), "\n")
