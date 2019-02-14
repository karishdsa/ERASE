#script to generate plots of the pvals at all alphas
#rm(list=ls())
library(ggplot2)
library(data.table)
BASEDIR <- "/home/kdsa/bootstrap/"

#source(paste0(BASEDIR, "scripts/randInteMultiRefs_functions.R"))
source(paste0(BASEDIR, "scripts/plots_functions.R"))

#read the pvals
pval_basedir <- paste0(BASEDIR, "results/")
pvals <- fetchAllPvalsAlphas (pval_basedir)
#saving pvals across tiss
acr_tiss <- pvals[ grep("ACROSS", pvals$ref), ]
outfile <- paste0(BASEDIR, "results/figures_data/data_acrossTissPvalsAllAlpha.tsv")
fwrite(acr_tiss, file=outfile , sep="\t")
rm(outfile)
#generating plots for each ref1 i.e. GWAS


transpvals <- fetchAllPvalsTrans(pval_basedir)
refList <- unique(transpvals$pvalsource)
refList <- refList[refList != "SPIDEX"]
outdir <- paste0(BASEDIR, "results/alpaPvalPlots/")
setwd(outdir)
ret <-lapply(refList, function(r){
        #r <- refList[1]
        outFile <- paste0(r,"_Alphapval.png")
        alp<- pvals[ grepl(r, ref),]
        tra <- transpvals[ grepl(r, ref),]
        return (plotEmpPvalAlphaWithIndiRefPval ( alp, tra, outFile))
        
      })

#individual plots for <tiss><ref1><ref2>
refList <- unique(pvals$ref)
ret <-lapply(refList, function(r){
  #r <- refList[1]
  outFile <- paste0(r,"_Alphapval.png")
  alp<- pvals[ grepl(r, ref),]
  #tra <- transpvals[ grepl(r, ref),]
  #return (plotEmpPvalAlphaWithIndiRefPval ( alp, tra, outFile))
  return (plotEmpPvalAlpha ( alp,  paste0("reverseX_",outFile))) #plotting reverese x axis and only bonferr pval
})
  
plotEmpPvalAlpha <- function( dtAlphaPvals, outFile){
  #function generates the plot of the alpha (in reverse 1 to 0) vs pvals in the datatable dt
  # for all the analysesin the ref column
  #no individual references pvals
  #dtAlphaPvals : with columns- alpha, pval, ref
    # plot saved as png in the current dir
  #dtAlphaPvals <- alp ; outFile <- paste0("reverseX_",outFile)
  
  dtAlphaPvals[ , neglog10pval:=(-log10(pval))]
  head(dtAlphaPvals, 2)
  
  
  bonpval <- 0.05/length(unique(dtAlphaPvals$alpha))
  title <- paste0("P-values at all alphas\n Bonferroni cut-off (orange line)= 0.05/",length(unique(dtAlphaPvals$alpha)), " = ", formatC(bonpval, format = "e", digits = 2) )
  #pval_pt5 <- 0.05
  #title <- paste0(title, "\n pval (violet line)= ", pval_pt5 )
  p <- ggplot( dtAlphaPvals, aes(x = alpha, y = neglog10pval) ) + 
    geom_point(size= 0.5 ,  shape=1, aes(color=ref) ) + 
    #geom_line() +
    geom_smooth(method = "loess", se = F, span=0.8,  lty ="dashed", size=0.25) +
    scale_x_reverse() +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  if (length(unique(dtAlphaPvals$ref)) ==1) {   #if only 1 analysis
    title <- paste( unique(dtAlphaPvals$ref), "-", title ) 
    #add the indi ref pval
    #+
    # geom_text(data=dtSingleRefPvals, aes(0, neglog10pval,label =pvalsource, vjust = -0.75), size=1) 
  } else {
    p <- p +
      facet_grid(ref~.) +
      theme(strip.text.y  = element_text(face="bold", size=2.5), 
            strip.background = element_rect(size=0.5)) +
      theme( axis.text  = element_text( size=6))
    
  } 
  
  p <- p +
    geom_hline(yintercept = -log10(bonpval), color="darkorange3" , linetype="dotdash", size=0.15) +
    geom_text(x=0, y= -log10(bonpval) + 0.15, size=2.5, aes(label = "Bonferroni cut off"),  vjust="inward",hjust="inward") +
    labs ( x= "Relative weight of GWAS",
           y = "-log10(pval)",
           title=title ) +
    theme(plot.title = element_text(size=SIZE_TITLE_TEXT-4, face="bold", hjust=0.5), 
          axis.title = element_text(size=6) ,
          legend.text=element_text(size=6),
          legend.title =element_text(size=6)) 
  
  ggsave(outFile,   plot = p)
  cat("\n Plot saved to ", outFile)
  return(T)
}



plotEmpPvalAlphaWithIndiRefPval <- function( dtAlphaPvals, dtSingleRefPvals, outFile){
  #function generates the plot of the alpha vs pvals in the datatable dt
  # for all the analysesin the ref column
  #also plots the individual references pvals for comparison
  #dtAlphaPvals : with columns- alpha, pval, ref
  #dtSingleRefPvals
  # plot saved as png in the current dir
  #dtAlphaPvals <- alp
  #dtSingleRefPvals <- tra
  
  dtAlphaPvals[ , neglog10pval:=(-log10(pval))]
  head(dtAlphaPvals, 2)
  dtSingleRefPvals[ , neglog10pval:=(-log10(pval))]
  
  
  bonpval <- 0.05/length(unique(dtAlphaPvals$alpha))
  title <- paste0("P-values at all alphas\n Bonferroni cut-off (orange line)= 0.05/",length(unique(dtAlphaPvals$alpha)), " = ", formatC(bonpval, format = "e", digits = 2) )
  pval_pt5 <- 0.05
  title <- paste0(title, "\n pval (violet line)= ", pval_pt5 )
  p <- ggplot( dtAlphaPvals, aes(x = alpha, y = neglog10pval) ) + 
    geom_point(size= 0.5 ,  shape=1, aes(color=ref) ) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  if (length(unique(dtAlphaPvals$ref)) ==1) {   #if only 1 analysis
    title <- paste( unique(dtAlphaPvals$ref), "-", title ) 
    #add the indi ref pval
    p <- p + geom_hline(data = dtSingleRefPvals, mapping=aes(yintercept = neglog10pval,color=pvalsource), 
                        linetype="dashed") #+
     # geom_text(data=dtSingleRefPvals, aes(0, neglog10pval,label =pvalsource, vjust = -0.75), size=1) 
  } else {
    p <- p +
      facet_grid(ref~.) +
      theme(strip.text.y  = element_text(face="bold", size=2.5), 
            strip.background = element_rect(size=0.5)) +
      theme( axis.text  = element_text( size=6))
    #add the indi ref pval
    p <- p + 
      geom_hline(data = dtSingleRefPvals, mapping=aes(yintercept = neglog10pval), 
                 linetype="dashed", size=0.25 ) +
      geom_text(data=dtSingleRefPvals, aes(0, neglog10pval,label =pvalsource, vjust = -0.75), size=1) 
  } 
  
  p <- p +
    geom_hline(yintercept = -log10(bonpval), color="darkorange3" , linetype="solid", size=0.25) +
    geom_hline(yintercept = -log10(pval_pt5), color="darkviolet" , linetype="solid", size=0.25) +
    labs ( x= "alpha",
           y = "-log10(pval)",
           title=title ) +
    theme(plot.title = element_text(size=SIZE_TITLE_TEXT-4, face="bold", hjust=0.5), 
          axis.title = element_text(size=6) ,
          legend.text=element_text(size=6),
          legend.title =element_text(size=6)) 
    
    #theme(plot.title = element_text(size=SIZE_TITLE_TEXT-4, face="bold", hjust=0.5), 
    #      axis.title = element_text(size=6) ,
    #      legend.text=element_text(size=6),
    #      legend.title =element_text(size=6)) +
    
  ggsave(outFile,   plot = p)
  cat("\n Plot saved to ", outFile)
  return(T)
}


fetchAllPvalsAlphas <- function(basedirPval){
  #returns all the alpha pvals in the basedirPval
  #can be a alphas for 1 or more analyses
  
  fPattern <- "AlphasPvals.tsv$"
  randFiles <- list.files(path=basedirPval, pattern = fPattern, recursive = T, full.names = T)
  recs <- lapply(randFiles, function(f){
    return( fread(f) )
  })
  recs <- rbindlist(recs)
  return(recs)
}


fetchAllPvalsTrans <- function(basedirPval){
  #returns all the pvals calculated using the Transformed scores for individual refs in the basedirPval
  
  fPattern <- "TransScoresPvals.tsv$"
  randFiles <- list.files(path=basedirPval, pattern = fPattern, recursive = T, full.names = T)
  recs <- lapply(randFiles, function(f){
    return( fread(f) )
  })
  recs <- rbindlist(recs)
  return(recs)
}

