###
#Variables
#----------
library(data.table)
library(ggplot2)
library(gridExtra)

#Variables for the plots
MULTI <- 1

TITLE_SIZE <- 10*MULTI
LEGEND_TITLE_SIZE <- 10*MULTI
LEGEND_TEXT_SIZE <- 8*MULTI
AXIS_TITLE_SIZE <- 10*MULTI
AXIS_TEXT_SIZE <- 8*MULTI
POSITION_DODGE_WIDTH <- 0.5

theme_clearBg <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"))

theme_legend <- theme(legend.title = element_text(colour="black", size=LEGEND_TITLE_SIZE, face="bold"), 
                      legend.text = element_text(colour="black", size=LEGEND_TEXT_SIZE))

theme_title <- theme(plot.title = element_text(size=TITLE_SIZE, face="bold", hjust = 0.5)) #hjust =0.5 to centre
theme_subtitle <- theme(plot.subtitle=element_text(size=TITLE_SIZE-1, hjust=0.5, face="italic"))
#X and Y axes labels, text
theme_axes <- theme(axis.title = element_text(face="bold", size=AXIS_TITLE_SIZE),
                    axis.text  = element_text( vjust=0.5, size=AXIS_TEXT_SIZE) )

fill_legend <- scale_fill_manual(values=c("turquoise4", "yellow3"))
fill_legend_meanavgRds <- scale_fill_manual(values=c("seagreen4", "slateblue3"))

###
#FUNCTIONS
#---------

##Function

#plot showing the distribution of mean number of SNPs per field2check bin e.g. avg.reads
#of the randomly selected SNPs [bins correspond to avg read depth] 
#across the 10K randomisations and obs

#Plot 2 - mean proportion of unique SNPs per bin
#returns the name of the file generated

plotMeanSNPsPerField2Chk4DistrBin <- function( obs,ran , field2chk4distr, fname.prefix, bins=NULL ) {
  #obs <- q1data.overlap; ran<- ran.from.q2data.overlap
  #bins <- NULL
  setDT(obs); setDT(ran)
  
  
  #get the avg no. of SNPs per bin across perm
  totPerPermBin <- ran[ , .(nSNPs.perbin.perm=length(cmp.col)), by="perm.no,bin.id"]
  avgPerBin <- totPerPermBin[ , .(nSNPs=mean(nSNPs.perbin.perm)), by="bin.id"]
  avgPerBin$type <- "Mean.nSNPs.perBin.acrossRand"
  setorder(avgPerBin, bin.id)
  
  #get the mean avg Reads per bin - Added att he twindUK data stage
  #avgPerBin <- ran[ , .(avg.reads=mean(avg.reads)), by=bin.id] 
  #avgPerBin$type <- "Mean.avgReads.perBin"
  meanAvgReadsPerPermBin <- ran[ , .(meanAvgReads.perbin.perm=mean(avg.reads)), by="perm.no,bin.id"]
  meanAvgReadsPerBin <- meanAvgReadsPerPermBin[ , .(meanAvgReads=mean(meanAvgReads.perbin.perm)), by="bin.id"]
  meanAvgReadsPerBin$type <- "Mean.AvgReads.perBin.acrossRand"
  setorder(meanAvgReadsPerBin, bin.id)
  
  #end get mean avg reads
  
  #mean proportion of unique SNPs per bin
  uniqPerPermBin <- ran[ , .(nUniqSNPs.perbin.perm=length(unique(cmp.col))), by="perm.no,bin.id"]
  head(uniqPerPermBin)
  head(totPerPermBin )
  setkeyv(uniqPerPermBin, c("perm.no", "bin.id"))
  setkeyv(totPerPermBin, c("perm.no", "bin.id"))
  res <- totPerPermBin[uniqPerPermBin ]
  head(res, 3)
  #totPerPermBin[bin.id==8 & perm.no==1,]
  #uniqPerPermBin[bin.id==8 & perm.no==1,]
  res$prop.uniq.SNPs <- res$nUniqSNPs.perbin.perm /res$nSNPs.perbin.perm
  avgPropUniqPerBin <- res[ , .(prop.uniq.SNPs.mean=mean(prop.uniq.SNPs)), by="bin.id"]
  #head(avgPropUniqPerBin)
  #(mean (res[ bin.id==92, prop.uniq.SNPs]))
  #avgPropUniqPerBin[bin.id ==92,]
  rm(res)
  
  #obs
  val<- obs[ ,.SD, .SDcol=field2chk4distr]
  if( ! is.null(bins)) {
    brk <- bins
  } else {
    if (BIN.MODE == "AVGREADS"){
      brk <- getBreaksBin(val ) 
    }
  }
  h <- hist(val$avg.reads, breaks=brk, plot=F)
  obsPerBin <- h$counts
  obsPerBin <- data.table( nSNPs = obsPerBin)
  #obsPerBin <- data.frame( avg.reads = obsPerBin, stringsAsFactors=F)
  obsPerBin$bin.id <- as.numeric(rownames(obsPerBin))
  setcolorder(obsPerBin, c( "bin.id", "nSNPs" ))
  obsPerBin$type <- "obs.nSNPs.perBin"
  
  
  #to get the avg reads per bin in obs
  #obsAvgReadPerBin <- h$breaks
  #obsAvgReadPerBin <- data.table( avgReadbrks = obsAvgReadPerBin) # using meanAvgReads as for rbindlist later need matching col names
  #obsAvgReadPerBin$bin.id <- as.numeric(rownames(obsAvgReadPerBin))
  setDF(val)
  #assing the bin id for the avg reads in the observed
  brk.probs <- getProbabilities(val$avg.reads, bins )
  obsmeanAvgReadsPerBin <- assignBinIndex (val , brk.probs ) #assign the bin ids for the avg reads using the breaks vals 
  #get the avg per bin
  setDT(obsmeanAvgReadsPerBin)
  obsmeanAvgReadsPerBin <- obsmeanAvgReadsPerBin[ , .(meanAvgReads=mean(avg.reads)), by="bin.id"]
  setorder(obsmeanAvgReadsPerBin, bin.id)
  obsmeanAvgReadsPerBin$type <- "obs.avgReads.perBin"
  #Kolmogorov Smirnov test - avgreads
  pval_avgreads <- ks.test( obsmeanAvgReadsPerBin$meanAvgReads, meanAvgReadsPerBin$meanAvgReads)$p.value
  
  
  #t_obs <- obsPerBin[bin.id >3, ]
  #Kolmogorov Smirnov test
  obsPerBin <- obsPerBin[ bin.id %in% obsmeanAvgReadsPerBin$bin.id, ] #added as the obserBin has the initial bins with 0 SNPs and affects the ks.test 
  pval <- ks.test( obsPerBin$nSNPs, avgPerBin$nSNPs)$p.value
  #pval <- 0.99
  
  comb <- rbindlist( list(avgPerBin, obsPerBin))
  #fwrite(comb, file="avgDistr.tsv", quote=F, sep="\t")
  stitle <- basename(fname.prefix)
  title <- paste0( "Distribution of mean #SNPs per bin (", field2chk4distr, ") of \n the randomly selected and the observed")
  title <- paste0(title, "\nKolmogorov Smirnov test p-value=", pval )
  #comb[ bin.id==92]
  
  p1 <-ggplot(comb, aes(x=bin.id, y=nSNPs, fill=type)) +
          ggtitle(title, subtitle = stitle) +
          geom_bar(stat="identity", position= position_dodge(width=POSITION_DODGE_WIDTH) ) +
          xlab("Bin ID") +
          ylab("#SNPs") +
          theme_axes +
          fill_legend +
          theme_legend +
          theme_title +
    theme_subtitle
  
  #plot for prop of unique SNPs per bin
  title <- paste0( "Distribution of mean proportion of unique SNPs per bin (", field2chk4distr, ") of \n the randomly selected set.")
  nrow( avgPropUniqPerBin)
  avgPropUniqPerBin$bin.id
  p2 <- ggplot(avgPropUniqPerBin, aes(x=bin.id, y=prop.uniq.SNPs.mean)) +
        geom_bar(stat="identity", fill="royalblue4", width=0.70) +
        ggtitle(title, subtitle = stitle) +
        xlab("Bin ID") +
        ylab("Mean proportion of unique SNPs per bin") +
        ylim(0,1) +
        theme_axes +
        theme_title +
    theme_subtitle
  #grid.arrange(p1, p2, nrow=1, ncol=1)       
  
  #plot for Mean avg reads per bin
  rm(comb)
  comb <- rbindlist( list(meanAvgReadsPerBin, obsmeanAvgReadsPerBin))
  title <- paste0( "Distribution of mean AvgReads per bin of \n the randomly selected and the observed")
  title <- paste0(title, "\nKolmogorov Smirnov test p-value=", pval_avgreads )
  p3 <-ggplot(comb, aes(x=bin.id, y=meanAvgReads, fill=type)) +
    ggtitle(title, subtitle = stitle) +
    geom_bar(stat="identity", position= position_dodge(width=POSITION_DODGE_WIDTH) ) +
    xlab("Bin ID") +
    ylab("Mean AvgReads") +
    theme_axes +
    fill_legend_meanavgRds +
    theme_legend +
    theme_title +
    theme_subtitle
  
  
  plots <- list(p1, p3, p2)
  
  outfile <- paste0(fname.prefix, "_DistrMeanSNPsPerBin.pdf")
  #ggsave(plot = p, filename = outfile, scale=2)
  ggsave( filename = outfile, marrangeGrob(grobs = plots, nrow=1, ncol=1), width=11, height=8.5)

#source("~/ase57/scripts/gg_multiplot.R")

#pdf(outfile, onefile=T, width=12, height=8 )
#par( mar=c(5,3,2,2)+0.1))
#multiplot( plots, cols=1, rows=1 )

#dev.off()

  return(outfile)
}

generatePlot_MeanSNPsPerField2Chk3DistrBin <- function(rdafile){
  #rdafile <- flist[2]
  load(rdafile)
  field2chk4distr <- "avg.reads"
  fname.prefix <- gsub( "_obsRandomVals.rda$", "", rdafile)
  return( plotMeanSNPsPerField2Chk4DistrBin ( q1data.overlap, ran.from.q2data.overlap , field2chk4distr, fname.prefix, bins=NULL )  )
  
}
## End function

