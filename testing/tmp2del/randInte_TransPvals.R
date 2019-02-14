#script : randInte_TransPvals.R
#Calculates the pvals for the ref1, ref2 using the transformed scores
#
#R --vanilla --slave < randInte_TransPvals.R tiss ref1 ref2
#                  tiss - the tissue ( e.g PUTM )
#                  ref1 - the GWAS dataset e.g SCZ
#                  ref2 - the 2nd ref dataset e.g SPIDEX, Annotation

#  Output generated in : <BASEDIR>results/<tiss>_<ref1>_<ref2>

tiss <- as.character( commandArgs()[4] )
ref1 <- as.character( commandArgs()[5] )
ref2 <- as.character( commandArgs()[6] )
#rm(list=ls())

# ref1 <- "SCZ2018";  ref2 <- "SPIDEX"; tiss <- "PUTM"

BASEDIR <- "/home/kdsa/bootstrap/"

source(paste0(BASEDIR, "scripts/randInteMultiRefs_functions.R"))
cat("\nStart ", date(), "\n")

inp_outdir <- paste0(BASEDIR, "results/", tiss, "_", ref1, "_", ref2)
setwd( inp_outdir)
refs <- paste0(ref1, "_", ref2)

randFilesPattern <- "*obsRandomValsMean_zscoreNmlDistr.rda"
randFiles <- list.files(pattern = randFilesPattern)


cat("\nTissue :", tiss, "\n")
cat("\nRef1 :", ref1, "\n")
cat("\nRef2 :", ref2, "\n")
nExpectedFiles <- 2

#check if 2 files found are for the 2 refs       
if ( (length(randFiles) != nExpectedFiles) | ( length(which( grepl (refs, randFiles))) != nExpectedFiles))  {
  cat("\n", nExpectedFiles," Transformed zscore files expected for 2 references.", length(randFiles)," files found\n")
} else {
  #calculate pval 
  trans_pval <- setDT(calculatePval4TransScores( randFiles))
  trans_pval[ , pvalsource := ifelse( grepl("GWAS",rda), ref1, ref2)]
  trans_pval$ref <- paste0(tiss, "_", refs)
  trans_pval$rda <- NULL
  outfile <- paste0(tiss, "_", ref1, "_", ref2, "_RandomisationEachRefTransScoresPvals.tsv")
  write.table(trans_pval, file=outfile,  sep="\t", quote=F, row.names=F)
  #save(alpha_pval,  file=outfile, compress=TRUE)
  cat( "\n\npvals saved to : ", outfile , " in ", getwd())
}

cat("\nEnd ", date(), "\n")
