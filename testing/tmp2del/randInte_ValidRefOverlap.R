# randInte_ValidRefOverlap.R
#Script gets the overlaps of the valid hets for the Tissue (PUTM, ACROSS_TISS (min FDR across tissue) ) with the references
#GWAS (ref1) SCZ, PD and
#ref2 - SPIDEX, Annotation - transcript specific

#Also generates plot for the -log10(pval) vs the zscores for the valid hets that overlaps with the GWAS summ data

#qrsh -l h_rt=4:0:0 -l tmem=7.9G,h_vmem=7.9G
#/share/apps/R/bin/R
#source("/home/skgtkds/randomisation/scripts/randInteMultiRefs_functions.R")
#BASEDIR <- "/SAN/neuroscience/ukbecvc/"
#outdir <- paste0(BASEDIR, "randomisation/files/overlapFiles/")

BASEDIR <- "/home/kdsa/bootstrap/"
source(paste0(BASEDIR, "scripts/randInteMultiRefs_functions.R"))
outdir <- paste0(BASEDIR, "files/overlapFiles/")
setwd(outdir)

tiss <- "PUTM"
#tiss <- "ACROSS_TISS"

#Get the avg read distributions

#avgreads <- read.delim( paste0(BASEDIR, "ase57/results/averageReadsPerHet_validSamples.tsv"), as.is=T, stringsAsFactors=F)  #252742
avgreads <- read.delim( paste0(BASEDIR, "files/ase_putmsnig/averageReadsPerHet_validSamples.tsv"), as.is=T, stringsAsFactors=F)  #252742
all.tissues <- c("PUTM", "SNIG")
#getting tiss specific cols and hetSNp
if(tiss %in% all.tissues ) {
  avgreads <- avgreads[ , names(avgreads)[ grep( paste0("hetSNP|", tiss), names(avgreads) )]] #252742
} else {
  avgreads <- avgreads[ , names(avgreads)[ grep( paste0(all.tissues, collapse="|"), names(avgreads), invert=T )]] #252742
}

names(avgreads) <- gsub( names(avgreads) [grep("Avg.Reads", names(avgreads) )] , "avg.reads", names(avgreads) )
names(avgreads) <- gsub( names(avgreads) [grep("min.FDR", names(avgreads) )] , "min.FDR", names(avgreads) )

avgreads <- avgreads[ which (!is.na(avgreads$min.FDR)), ]#196049 - PUTM, 252742 - across
avgreads$cmp.col <- avgreads$hetSNP

cat("\n", tiss, "\n valid hets=",  nrow(avgreads) )
cat( "\n ASEs =",  length(which(avgreads$min.FDR < 0.05)) )
cat( "\n nonASEs =",  length(which( !is.na(avgreads$min.FDR) & avgreads$min.FDR >= 0.05)), "\n\n" )

#hetrsidfile <- paste0( BASEDIR, "ase57/results/gwas_studies/validHets_rsid.tsv")
hetrsidfile <- paste0( BASEDIR, "files/ase_putmsnig/gwas_studies/validHets_rsid.tsv")

#Get the overlap for the different references

#ref1 data- GWAS summary datasets
#W2H2015
refFile <-  refFile <- "/home/kdsa/bootstrap/files/W2H2015/GIANT_2015_WHRadjBMI_COMBINED_EUR_hg19.tsv"
refdata <-"W2H2015"
gwas <- getGWASrecs (refFile)
ref1 <- gwas
rm(gwas)


#RA2014
refFile <- "/data/LDScore/GWAS/RA2014/RA_GWASmeta_European_v2.txt.gz"
refdata <-"RA2014"
gwas <- getGWASrecs (refFile)
colnames(gwas)[colnames(gwas)=="Position(hg19)"] <- "Position"
gwas$cmp.col <- paste0(gwas$Chr, ":", gwas$Position )
gwas <- updatepvaleq0 (gwas)
ref1 <- gwas
rm(gwas)

#SCZ
scz <- "/SAN/neuroscience/WT_BRAINEAC/GWAS/scz/sczfullsnpsresults.txt.gz"
refdata <-"SCZ"
gwas <- getGWASrecs (scz)
gwas$cmp.col <- paste0( gsub("chr", "", gwas$hg19chrc ), ":", gwas$bp)
ref1 <- gwas
rm(gwas)

scz <- "/SAN/neuroscience/ukbecvc/ase57/files/SCZ2018/clozuk_pgc2.meta.sumstats.txt.gz"
refdata <- "SCZ2018"
delim <- "space"
gwas <- getGWASrecs (scz, delim)
gwas$cmp.col <- paste0(gwas$CHR, ":", gwas$BP )
gwas$p <- gwas$P
gwas <- updatepvaleq0 (gwas)
ref1 <- gwas
rm(gwas)

scz <- "/SAN/neuroscience/ukbecvc/ase57/files/SCZ2013/scz.swe.pgc1.results.v3.txt.gz"
refdata <- "SCZ2013"
ref1 <- getGWASrecs (scz) #default = tab
ref1$cmp.col <- paste0(ref1$hg19chr, ":", ref1$bp )
ref1 <- updatepvaleq0 (ref1)
#ref1 <- gwas
#rm(gwas)


scz <- "/SAN/neuroscience/ukbecvc/ase57/files/SCZ2011/pgc.scz.full.2012-04.tsv"
refdata <- "SCZ2011"
ref1 <- getGWASrecs (scz) #default = tab  #1252901 recs hg18, all have rsids
ref1$cmp.col <- ref1$snpid
ref1 <- updatepvaleq0 (ref1)

hetrsid <- read.delim(hetrsidfile, as.is=T, stringsAsFactors=F)
hetrsid <- hetrsid[ , c("hetSNP", "rsid")]#244088
#add the avg reads and the minFDR
avgreads <- merge(hetrsid, avgreads, by="hetSNP")
avgreads$cmp.col <- avgreads$rsid


#PD
#pd <- "/SAN/neuroscience/ukbec/common_files/GWAS/PD_13_03_17/AllResults.txt.zip"  ##7893273
#refdata <-"PD"

pd <- "/SAN/neuroscience/ukbecvc/ase57/files/PD2017_meta5/resultsForCojo_april17th2017.tab.gz"  #8055803
refdata <-"PD2017_Meta5"

pd <- "/SAN/neuroscience/WT_BRAINEAC/GWAS/pd/pdMeta1-forLDSC.txt"
refdata <-"PD2014"

pd <- "/SAN/neuroscience/ukbecvc/ase57/files/PD2017_Chang_ex23andMe/META_no231_formatted.txt"
refdata <-"PD2017ChangEx23andMe"   #10756947 recs

gwas <- getGWASrecs (pd)  #


if(refdata=="PD") {
  gwas$cmp.col <- paste0( gsub("chr", "", gwas$Chr ), ":", gwas$Bp)
  #names(gwas)[ grep("P.value.1", names(gwas) ) ] <- "p"
  #names(gwas)[ grep("P.value", names(gwas) ) ] <- "p"
  gwas$p <- gwas$P.value
  gwas <- updatepvaleq0 (gwas)  #is called in getGWASrecs- however due to issue with 2 P.value fields, p value field is being set here.
  # Therefore calling function here again.
} else {
  if ( refdata== "PD2017_Meta5" ) {
    gwas$cmp.col <- gwas$MarkerName
  } else {
    if ( refdata== "PD2014" ) {
      gwas$cmp.col <- gsub("^chr", "", gwas$MarkerName ) # also contains recs with rsid
      #get valid hets with rsid and add to avg reads
      hetrsid <- read.delim(hetrsidfile, as.is=T, stringsAsFactors=F)
      hetrsid <- hetrsid[ , c("hetSNP", "rsid")]#244088
      #add the avg reads and the minFDR
      hetrsid <- merge(hetrsid, avgreads, by="hetSNP", all.x=T)
      hetrsid$cmp.col <- hetrsid$rsid
      
    } else {
      if ( refdata == "PD2017ChangEx23andMe") {
        gwas$cmp.col <- gsub("chr", "", gwas$MarkerName ) # no recs with rsid
      }
    }
  }
}


ref1 <- gwas
rm(gwas)

#AD
ad <- "/SAN/neuroscience/WT_BRAINEAC/GWAS/ad/IGAP_stage_1.txt"
refdata <- "AD"
gwas <- getGWASrecs (ad)

gwas$cmp.col <- paste0( gwas$Chromosome, ":", gwas$Position)
ref1 <- gwas
rm(gwas)



#common
#valid overlap with ref
ref1.field4mean <- "p"
avg.ref1 <- getHetOverlapRef(avgreads, ref1, ref1.field4mean )  #PUTM 174470 overlap ref1SCZ, 165733 PD
if ( refdata== "PD2014" ) {
  #get the overlap of the rsids
  rsid.ref1 <- getHetOverlapRef(hetrsid, ref1, ref1.field4mean )
  rsid.ref1$rsid <- NULL
  comb <- rbind(rsid.ref1, avg.ref1)
  avg.ref1 <- comb
  rm(comb)
  avg.ref1$cmp.col <- avg.ref1$hetSNP
}

minrecs <- getMinPval ( avg.ref1)  # SNPs with > 1 pval, min pval for the SNP is assigned
avg.ref1$p <- NULL
avg.ref1 <- unique ( avg.ref1 )        #across - 227751 overlap SCZ; 217219 overlap PD; 224126 Meta5
avg.ref1 <- merge( avg.ref1, minrecs, by ="cmp.col")
avg.ref1$ref <- refdata
valid.gwas <- avg.ref1
rm(avg.ref1)

#saving the valid hets that overlap with the ref
fname <- paste0(tiss, "_validhets_", refdata, "_Overlap.rda" ) #PUTM_validhets_SCZ_Overlap.rda
save(valid.gwas, file=fname, compress=TRUE)


cat("\n No. of Overlapping SNPs=", nrow(valid.gwas))
cat("\n No. ASEs", length( which(valid.gwas$min.FDR < 0.05)) )
cat("\n No. non ASEs", length( which(valid.gwas$min.FDR >= 0.05)) )

#checkign the disctribuion of Avg reads in the ASEs vs non ASEs
source(paste0(BASEDIR,"scripts/process_ASEdatasets_functions.R") )
bwidthlist <- 2
generatePlot_AvgReadDistribution4Bin(fname,bwidthlist)

##################################

#get the overlap of the 2 refs and valid hets ( intersection )

setwd("/SAN/neuroscience/ukbecvc/ase57/results/randomisationInte/overlapFiles/")
#scz, spidex
ref1.rda <- "ACROSS_TISS_validhets_SCZ_Overlap.rda"
ref2.rda <- "ACROSS_TISS_validhets_SPIDEX_Overlap.rda"
out.rda <- "ACROSS_TISS_validhets_SCZ_SPIDEX_Overlap.rda"

ref1.rda <- "ACROSS_TISS_validhets_SCZ2018_Overlap.rda"
ref2.rda <- "ACROSS_TISS_validhets_SPIDEX_Overlap.rda"
out.rda <- "ACROSS_TISS_validhets_SCZ2018_SPIDEX_Overlap.rda"

#PD, spidex
ref1.rda <- "ACROSS_TISS_validhets_PD_Overlap.rda"
ref2.rda <- "ACROSS_TISS_validhets_SPIDEX_Overlap.rda"
out.rda <- "ACROSS_TISS_validhets_PD_SPIDEX_Overlap.rda"

ref1.rda <- "ACROSS_TISS_validhets_PD2017ChangEx23andMe_Overlap.rda"
ref2.rda <- "ACROSS_TISS_validhets_SPIDEX_Overlap.rda"
out.rda <- "ACROSS_TISS_validhets_PD2017ChangEx23andMe_SPIDEX_Overlap.rda"


ref1.rda <- "ACROSS_TISS_validhets_PD2017_Meta5_Overlap.rda"
ref2.rda <- "ACROSS_TISS_validhets_SPIDEX_Overlap.rda"
out.rda <- "ACROSS_TISS_validhets_PD2017_Meta5_SPIDEX_Overlap.rda"

#common
#calling function to get the intersection - for diff ref1s
getIntersectionValRefs(ref1.rda, ref2.rda, out.rda )


###################################

#Generating plot for the -log10(pval) vs the zscores for the valid hets that overlaps with the GWAS summ data

#SCZ
f.rda <- "PUTM_validhets_SCZ_Overlap.rda"

#PD
f.rda <- "PUTM_validhets_PD_Overlap.rda"

load(f.rda)
p.fname <- gsub(".rda", "_NeglogpvalVsZscorePlots.pdf", f.rda)

generatePlot_ZvsNeglogpval( valid.gwas, p.fname)


#ref2 data
#spidex
#get the overlap with valid hets

refdata="SPIDEX"
spi <- getrecsSpiValOverlap()  #42979
spi <- spi[ , c("hetSNP", "dpsi_max_tissue", "dpsi_zscore")]
spi$abs_dpsi_max_tissue <- abs( spi$dpsi_max_tissue)

spi <- merge(spi, avgreads, by = "hetSNP")  #37,000 - PUTM ; 42,979 - across tiss

spi$cmp.col <- spi$hetSNP
spi$ref <- refdata

valid.ref <- spi
rm(spi)

#annotation- transcript specificity

#saving the valid hets that overlap with the ref
fname <- paste0(tiss, "_validhets_", refdata, "_Overlap.rda" ) #PUTM_velidhets_SPIDEX_Overlap.rda ; ACROSS_TISS_validhets_SPIDEX_Overlap.rda
save(valid.ref, file=fname, compress=TRUE)


#############Functions######################
source("/home/skgtkds/ase57/scripts/gwasDiseaseDatasets_controlSize_functions.R")
generatePlot_ZvsNeglogpval <- function( df, p.fname) {
  df.z <- convertPval2Z(df)
  df.z$neglogp <- -log10(df.z$p)
  
  pdf(file=p.fname)
  par(mfrow = c(2,2))
  title <- gsub(".pdf", "", p.fname)
  title.size <- 0.75
  
  
  
  ylabel <- "Z score:(sqrt(qchisq(pval, 1, lower.tail=F) ) )"
  maxval <- max(df.z$Z, df.z$neglogp )
  lim <- c(0, maxval+1)
  plot( df.z$neglogp, df.z$Z , xlab="-log10(p-val)", ylab=ylabel, xlim=lim, ylim=lim, main=paste0( title, "\n validhets" ), cex.main=title.size  )
  abline(0,1, col="red", lty=2 )
  #ases
  ase <- df.z[ which(df.z$min.FDR < 0.05), ]
  maxval <- max(ase$Z, ase$neglogp )
  lim <- c(0, maxval+1)
  plot( ase$neglogp, ase$Z , xlab="-log10(p-val)", ylab=ylabel,  xlim=lim, ylim=lim, main=paste0( title, "\nASEs" ), cex.main=title.size )
  abline(0,1, col="red", lty=2 )
  #nonASEs
  
  ase <- df.z[ which(df.z$min.FDR >= 0.05), ]
  maxval <- max(ase$Z, ase$neglogp )
  lim <- c(0, maxval+1)
  plot( ase$neglogp, ase$Z , xlab="-log10(p-val)", ylab=ylabel, xlim=lim, ylim=lim, main=paste0( title, "\nNon ASEs" ), cex.main=title.size )
  abline(0,1, col="red", lty=2 )
  dev.off()
  cat("\nOutput file: ", p.fname)
  return(T)
}

