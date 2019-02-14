#functions, variables for plots generated
library(ggplot2)
#Functions
#---------
SIZE_GEOM_POINT <- 1.5
SIZE_X_AXIS_TEXT<-12
SIZE_TITLE_TEXT <-12
theme_title <- theme(plot.title = element_text(size=SIZE_TITLE_TEXT, face="bold", hjust=0.5))


#Mainly for simVaryAseNum_plots.R
point_col <- "blue"

fetchPvals <- function(srcDir) {
  #returns the p-vals for the simulations for the ref
  # if the number of pval files is not = to the number of simulations in SIMS, returns a message
  #srcDir <- dir  
  cat("\n Getting simulation p-vals from - ", srcDir, "\n")
  
  fnames <- list.files(path=srcDir, pattern="*obsRandomValsZscorePval.rda",full.names = T, recursive=T)
  head(fnames)
  
  # if ( length(fnames) != SIMS ) {
  #    msg <- paste0("No. files != ", SIMS)
  #    return ( msg) 
  #  } 
  
  #fnames <- head(fnames) 
  pvals <- lapply(fnames, function(filename){
    #filename <- fnames[1]
    load(filename)
    if(length( grep("^gwasPval", basename(filename))) ==0){
      #if gwasPval, then the lower tail T pnorm pval was calculated
      pnorm.val.lowerTail.T <- NA
    }
    return(data.table(filename,pnorm.pval, pnorm.val.lowerTail.T))
  })
  head(pvals)
  pvals <- rbindlist(pvals)
  gc()
  
  return(pvals)
}

#returns a data table with the nASes for the percentage list provided
getNase4Pct <- function(tot, lstpct){
  
  lstnASEs <- lapply(lstpct, function(pct){
    nASEs <- round ((pct/100)*tot )
    return(data.table(pct, nASEs))
  })
  dtnASEs <- rbindlist(lstnASEs)
  return(dtnASEs)
} 

#function generates the distribution of the pvals obtained with the varying no. of ASEs
#with the differnt column used to rank the SNps i.e. gwas pval and neglog10pval
plotVaryAseNum_varyrankBy <- function(dt, title, outFile){
  #  dt <- recs
  setorder(dt, pct)
  ggplot( dt, aes(x = pct, y = -log10(pnorm.pval)) ) +
    ggtitle(title) +
    geom_boxplot( aes(group = nAse), color = point_col, fill="white") +
    geom_jitter( size = SIZE_GEOM_POINT, color = point_col, 
                 position=position_jitter() ) +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12)) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")+
    scale_x_continuous(breaks=unique(dt$pct),
                       labels=unique(dt$xlabel)) +
    xlab("Percentage(nASEs)") +
    theme_title +
    facet_wrap(lowertail~rankBy, ncol=2,  labeller=label_both)
  #    facet_grid(lowertail~rankBy,  labeller=label_both) # displays blank plot for neglog10 lowertail T 
  
  ggsave(outFile,width = 20, height = 10)  
}

#function generates the distribution of the pvals obtained with the varying no. of ASEs
plotVaryAseNum <- function(dt, title, outFile){
  #outFile <- outfname
  #head(dt)
  setorder(dt, pct)
  yint <- 0.05/nrow(dt) #Bonferroni cutoff
  title <- paste0(title, "\n Red line: Bonferroni cutoff- 0.05/", nrow(dt))
  #yint <- 0.05
  ggplot( dt, aes(x = pct, y = -log10(pnorm.pval)) ) +
    ggtitle(title) +
    geom_boxplot( aes(group = nAse), color = point_col, fill="white") +
    geom_jitter( size = SIZE_GEOM_POINT, color = point_col, 
                 position=position_jitter() ) +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12)) +
    geom_hline(yintercept=-log10(yint), linetype="dashed", color = "red")+
    scale_x_continuous(breaks=unique(dt$pct),
                       labels=unique(dt$xlabel)) +
    xlab("Percentage(nASEs)") +
    theme_title
  
  ggsave(outFile, width = 10, height = 7 )
}

#Plot the FDRs
#function generates the distribution of the FDRs obtained with the varying no. of ASEs
#cut off at 5% and 1% FDR
plotVaryAseNum_FDR <- function(dt, title, outFile){
  #outFile <- outfname
  #head(dt)
  setorder(dt, pct)
  
  yint <- 0.05
  title <- paste0(title, "\n Red line: FDR 5%, Black line = FDR 1%")
  ggplot( dt, aes(x = pct, y = -log10(FDR)) ) +
    ggtitle(title) +
    geom_boxplot( aes(group = nAse), color = point_col, fill="white") +
    geom_jitter( size = SIZE_GEOM_POINT, color = point_col, 
                 position=position_jitter() ) +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12)) +
    geom_hline(yintercept=-log10(yint), linetype="dashed", color = "red")+
    geom_hline(yintercept=-log10(0.01), linetype="dashed", color = "black")+
    scale_x_continuous(breaks=unique(dt$pct),
                       labels=unique(dt$xlabel)) +
    xlab("Percentage(nASEs)") +
    theme_title
  
  ggsave(outFile, width = 10, height = 7 )
}

#For RandomisationRefeQTLs_plots
SIZE_LINE_YINT<-0.25
SIZE_GEOM_LINE<-0.50

#-------End randomisationRefeQTLs_plots

#Functions for simulations-Type 1 error
#-------------------------------------
fetchSimPvals <- function(srcDir, nSims) {
  #returns the p-vals for the simulations for the ref
  # if the number of pval files is not = to the number of simulations in SIMS, returns a message
  
  cat("\n Getting simulation p-vals from - ", srcDir, "\n")
  
  fnames <- list.files(path=srcDir, pattern="*obsRandomValsZscorePval.rda",full.names = T)
  if ( length(fnames) != nSims ) {
    msg <- paste0("No. files != ", nSims)
    return ( msg) 
  } 
  
  #fnames <- head(fnames) 
  pvals <- lapply(fnames, function(f){
    #f <- fnames[1]
    load(f)
    dt <- data.table(pnorm.pval)
    dt$filename <- f
    return(dt)
  })
  gc()
  #return(unlist(pvals))
  return(rbindlist(pvals))
}

getSimPvals4AllSrcInDir <- function(maindir) {
  #get the pvals for all the sources for a set of simulations 
  #e.g. maindir= /home/kdsa/bootstrap/results/simulation_validhetsAvgRdChk  
  #returns a data.table with colnames, simpvals, source
  
  #get the sources the sim was run on - from the directory names
  dirlist <- dir(maindir, full.names = T) 
  
  pvals <- lapply( dirlist, function(simdir){
    #simdir <-dirlist[1]
    # fetch the data
    cat("\n Sim Dir :", simdir , "\n")
    r <- basename(simdir)
    cat("Source :", r , "\n")
    
    simpvals <- fetchSimPvals( simdir, nSims)
    #head(simpvals)
    dtpvals <- data.table( simpvals) 
    rm(simpvals)
    dtpvals <- dtpvals[ , source:=r]
    dtpvals$simdir <-simdir
    #head(dtpvals)
    
    cat("\n nSims where pvals < 0.05 =")
    cat( nrow (dtpvals[ (simpvals < 0.05)]), "\n" )
    return(dtpvals)
  } )
  pvals <- rbindlist(pvals)
  return(pvals )
}



#End functions
