#common functions used while processing the UKBEC, Twins UK and dermitzakis LCL data

library(ggplot2)
library(gridExtra)
#functions
#---------
getUKBEC_AvgReadsForTiss <- function(tiss) {
  #returns the avg reads, min.FDR for the tissue
  #tiss = PUTM, SNIG, ACROSS_TISS
  
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
  
  return(avgreads)
  
}

generateMinFDR_AvgReads_ASEfiles <- function(f, filter=T){
  
  #retval <- lapply(flist, function(f, filter=T) {
  #f <- flist[1] #"EB_F_ase"
  
  dt <- fread(f)
  head(dt, 2)
  cat("\nFile : ", f )
  cat( "\nNo. of samples = ", length(unique(dt$SampleID)) ) #716 Fat
  dt$cmp.col <- paste0(dt$CHR, ":", dt$POS)
  cat("\n No. of unique sites across samples = ", length(unique(dt$cmp.col))) #no. of sites per sample - before and after filtering
  
  if(filter==T){
    #filter at total reads =30 reads and both alleles observed
    dt<- dt[ TOTAL_COUNT >=30  & BOTH_ALLELES_SEEN ==1]
    cat("\n\nFiltering for total read count >=30 and both alleles seen \n")
    cat("\n No. of unique sites across samples = ", length(unique(dt$cmp.col))) #no. of sites per sample - before and after filtering
    
  }
  #calculate FDR per sample
  dt[ , FDR:= p.adjust(PVALUE, method="fdr"), by="SampleID"] 
  #head(dt$FDR) # 1.000000e+00 6.705123e-01 2.476051e-73
  #chk :   #sam <- dt[ SampleID=="EB53021" ,]
  #head( p.adjust(sam$PVALUE, method="fdr")) # [1] 1.000000e+00 6.705123e-01 2.476051e-73 
  
  #get the minFDR across samples
  dt[ , min.FDR:= min(FDR), by="cmp.col"] 
  head(dt, 2)   #min(dt[ cmp.col=="11:2150444", FDR])
  cat("nhetSNPs= ", length(unique(dt$cmp.col)), "\n") 
  cat("nSig min.FDR < 0.10= ", length( unique(  dt[ min.FDR <0.10, cmp.col])), "\n")
  cat("nSig min.FDR < 0.05= ", length( unique(  dt[ min.FDR <0.05, cmp.col])), "\n")
  
  #get Avg read depth across samples
  dt[ , avg.reads:= mean(TOTAL_COUNT), by="cmp.col"]
  head(dt, 2)   #mean( dt[ cmp.col=="11:2150444", TOTAL_COUNT])
  
  dt <- unique( dt[ , .(RSID, cmp.col, min.FDR, avg.reads)]) 
  outfile <- paste0(f, "_minFDR_avgreads.tsv")
  if(filter==F){
    outfile <- paste0("noFilter_", outfile)
  }
  fwrite(dt, file=outfile, quote=F, sep="\t")
  
  #})
}


plotAvgReads <- function(qrOverlap, title, outpdf, bwidth) {
  #plots the average reads  and zooms in
  # the data table qrOverlap column avg.reads is plotted in the hist
  # colour is by the column type i.e. ASE not ASE
  #pdf(file=outpdf)
  
  #plots
  p <- ggplot(qrOverlap, aes(x=avg.reads, fill=type)) +
    geom_histogram(binwidth=bwidth, position="dodge") + labs(x="average reads across samples at a hetSNP" ) +
    ggtitle(title)
  head(qrOverlap)
  
  #zooming
  dt <- qrOverlap[avg.reads<500,]
  ztitle <- paste0("Zooming -avgreads <500 \n", title)
  p1 <- ggplot(dt, aes(x=avg.reads, fill=type)) +
    geom_histogram(binwidth=bwidth, position="dodge") + labs(x="average reads across samples at a hetSNP" ) +
    ggtitle(ztitle)
  
  dt <- dt[avg.reads>50,]
  ztitle <- paste0("Zooming -avgreads >50 & <500 \n", title)
  p2 <- ggplot(dt, aes(x=avg.reads, fill=type)) +
    geom_histogram(binwidth=bwidth, position="dodge") + labs(x="average reads across samples at a hetSNP" ) +
    ggtitle(ztitle)
  
  plist <- list(p, p1, p2)
  ggsave(outpdf,   marrangeGrob(grobs = plist, nrow=1, ncol=1), width=11, height=8.5)
  #dev.off()
  
}

generatePlot_AvgReadDistribution4Bin <- function(f,bwidthlist, refSingle="T"){
  #f= file to be read containg col names avg.reads, min.FDR
  #refSingle = T - single regerence univariate method; "F" - multiple references - multivariate method
  cat("\n", f)
  
  #if overlap rda file
  if(grepl("Overlap.rda",f )) {
    load(f)
    if(refSingle=="T") {
      recs <- valid.gwas
    } else {
      recs <- valid.ref
    }
    setDT(recs)
  } else {
    recs <- fread(f)
  }
  
  head(recs,3)
  recs$type <- "Not ASE" 
  recs[ min.FDR < 0.05, type :="ASE"]
  print(  table(recs$type))
  
  ase <- recs[ type=="ASE", .(avg.reads)]
  head(ase$avg.reads)
  notase <- recs[ type=="Not ASE", .(avg.reads)]
  cat("\n Range of Avg reads : ASE=", range(ase$avg.reads), "; Not ASE =", range(notase$avg.reads))
  ks_pval <- ks.test(  ase$avg.reads, notase$avg.reads)$p.value
  
  title <- paste0(basename(f), "\n- Distribution of average reads across samples at a hetSNP\n Kolmogorov smirnov test p-val:",ks_pval )
  
  outprefix <- gsub(".rda|.tsv", "", basename(f))
  ret <- lapply(bwidthlist, function(x) {
    
    #outfile <- gsub(".rda", "_DistrAvgReads.pdf", basename(rdafile))
    
    outfile <- paste0(outprefix, "_bin", x, "_DistrAvgReads.pdf")
    plotAvgReads(recs, title, outfile,x)
    cat("\n", outfile)
  })
  
}


#
readLappOriFile <- function() {
  #function reads the  Lappalainen file containing the ASE data
  # returns a data.table with the added column cmp.col
  fname <- paste0(BASEDIR,"/files/dermitzakis/dl/GD462.ASE.COV8.ANNOTPLUS.SINFO.txt")
  dt <- fread(fname, fill = T)
  head(dt,2) # the the col names need to be shifted by 2 to the left.- has duplicate col names
  names(dt)[1] <- "rm1"
  names(dt)[2] <- "rm2"
  names(dt)[3] <- "SAMPLEID"
  setnames(dt, old = names(dt), new <- c( names(dt)[3: length(names(dt))], c( "tmp1", "tmp2")) )
  
  names(dt)
  head(dt,2)
  
  #check that the sampleid and sample col data are the same
  stopifnot(identical(dt$SAMPLE, dt$SAMPLEID))
  
  dt$tmp1 <- NULL
  dt$tmp2 <- NULL
  nrow(dt) #6284377
  dt$cmp.col <- paste0(dt$CHR, ":", dt$POS)
  length(unique(dt$cmp.col)) #352852
  length(unique(dt$SAMPLE)) #462
  return(dt)
}


getAvgMinMaxRatios <- function(dt) {
  #returns the avergae, min and max ratio across samples for a het
  #expects col names 'cmp.col' = group byu this column across samples e.g. hetSNP chr:pos
  #                   'ratio'
  #dt <- ratios
  vals <- dt[, .(avg.ratio = mean(ratio, na.rm = T), min.ratio = min(ratio, na.rm = T), max.ratio = max(ratio, na.rm = T)), by = cmp.col ]
  return(vals)
}

calculateAllelicRatios <- function( dt){
  #function orders the alleles alphabetically and calculates the ratios
  #(c1++0.5)/ (c2+0.5)
  
  #swap alleles that are not alphabetically sorted and the corresponding counts
  dt[ allele1>allele2,  c("allele1", "allele2", "count1", "count2") := .(allele2, allele1, count2, count1)]
  dt[ , ratio := ( count1 + 0.5) / ( count2 + 0.5  )]
  return(dt)
}

plotRatios <- function(dt, title, outpdf, plottype, bwidth=2) {
  #plots the allelic ratios (counts) and zooms in
  # the data table dt column avg.ratio  is plotted in the hist
  # colour is by the column type i.e. ASE not ASE
  #pdf(file=outpdf)
  #plot type = hist, density
  
  #dt <- recs;head(dt); bwidth <-2; outpdf <- outfile
  #plottype <- "density"
  
  #all data irrespective of ASE, non ASE
  m_title <- paste0("Distribution of allelic ratios \n All sites - Range of alleleic ratios ", 
                    paste0 ( formatC(range(dt$avg.ratio), format="e", digits=2), collapse=" - " ) )
  
  
  p_all <- ggplot(dt, aes(x=avg.ratio )) +
          geom_density() +
          ggtitle(m_title)
  
  #by ASEs and non ASEs
  rangeVals <-getRatioRange(dt)
  m_title <- paste0(title, "\nRange of alleleic ratios for ASEs : ", rangeVals[1], 
                  "; for non ASEs : ",rangeVals[2])
  
  cat("\n", m_title)
  if( plottype== "density") {
    p <- ggplot(dt, aes(x=avg.ratio, color=type)) +
      geom_density() 
    
  } else { #hist
    p <- ggplot(dt, aes(x=avg.ratio, fill=type)) +
      geom_histogram(binwidth=bwidth, position="dodge")  
      ##labs(x="average reads across samples at a hetSNP" ) +
    
    }
  
  p <- p +
    ggtitle(m_title)
  
  head(dt)
  #zooming
  cutoffLst <- c( 500, 100, 50, 10)
  pLst <- list()
  pLst <- lapply(cutoffLst, function(cutoff){
    #cutoff <- cutoffLst[5]
    z_dt <- dt[avg.ratio<cutoff,]
    rangeVals <-getRatioRange(z_dt)
    rangetitle <- paste0("\nRange of alleleic ratios for ASEs : ", rangeVals[1], 
                         "; for non ASEs : ",rangeVals[2])
    
    ztitle <- paste0("Zooming -avg ratio < ", cutoff, rangetitle, "\n", title)
    cat("\n avg.ratio <", cutoff, " ", rangetitle)
    if(plottype == "density") {
      return( ggplot(z_dt, aes(x=avg.ratio, color=type)) +
                geom_density() + 
                ggtitle(ztitle) )
    } else {
      return ( ggplot(z_dt, aes(x=avg.ratio, fill=type)) +
        geom_histogram(binwidth=bwidth, position="dodge") + 
        #labs(x="average reads across samples at a hetSNP" ) +
        ggtitle(ztitle) )
    }
  })
  plist <- c(list(p_all, p), pLst)
  #length(plist)
  ggsave(outpdf,   marrangeGrob(grobs = plist, nrow=1, ncol=1), width=11, height=8.5)
  cat("\n Plots saved to ", outpdf)
  return(outpdf) 
}

#to get the range of ratios in the ASEs and non ASEs
getRatioRange <- function(dt_ratio) {
  #returns the range as a list ASEs, non ASEs
  #rangeValsASE <- paste0( range( dt_ratio[ type=="ASE", .(avg.ratio)])  , collapse=" to ")
  #rangeValsNotASE <- paste0( range( dt_ratio[ type=="NotASE", .(avg.ratio)])  , collapse=" to ")
  rangeValsASE <-  paste0( formatC(range( dt_ratio[ type=="ASE", .(avg.ratio)]), format="e", digits=2) , collapse=" to ")
  rangeValsNotASE <- paste0( formatC( range( dt_ratio[ type=="NotASE", .(avg.ratio)]), format="e", digits=2) , collapse=" to ")
  
  return( list(rangeValsASE, rangeValsNotASE) )
}

