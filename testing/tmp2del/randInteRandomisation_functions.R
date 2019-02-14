#Script : randInteRandomisation_functions.R
# Functions to run randomisation


source(paste0(BASEDIR,"/scripts/randInte_Plotsfunctions.R"))

####### Functions for processing pre-randomisation ########
source(paste0(BASEDIR,"/scripts/randInte_ValidRefOverlap_functions.R"))

updatepvaleq0 <- function(df) {
  #Assign new min if pval==0
  if( length(which(df$p==0) ) >0 ) {
    next.min <- min( df[ which(df$p != 0),"p" ] )
    new.min <- next.min*1e-100
    df[ which(df$p==0), "p"] <- new.min
  }
  return(df)
}

getMinPval<- function(df) {
  #function returns a dataframe with the min pval per SNP, if a SNP has > 1 pval

  recs <-  data.frame( p =  with( df,( tapply( p, cmp.col, function(z) { min(z, na.rm=T ) })) ) , stringsAsFactors=F)
  recs$cmp.col <- rownames(recs)
  recs$p <- as.numeric( recs$p )
  rownames(recs) <- NULL
  return(recs)
}



#######Randomisation ##########
library(TeachingDemos)
RSID.REPEAT <-TRUE

#Internal setting for bins
BINWIDTH1 <- 2
BINWIDTH1_END <- 200
BIN.MODE <- "AVGREADS"  # if the bin list is provided the use that , else for now if AVGREADS (set as default) then
                        # use binwidth 2 for avgreads  <200 and remaining in 1 bin
                        # can add on if internally want to set for MAF instead of passing the bins

SEEDSTRING <- "RANDOMISATION"
#seed for the random number generator from a character string
# if no char string is passed, SEEDSTRING is used as default

getSeed4chrStr<- function(charval=SEEDSTRING){
  cat("\nSeed string = ", charval, "\n")
  return( char2seed(charval, set=F) )


}

# randomisation ref - call functions  if distribution check to be done on the randomly selected
#set of SNPs, such that thhey match that of the query 1 data (Category SNPs)
#q1data.overlap <- q1.ref ; q2data.overlap <- q2.ref; field4mean <- ref.field4mean ; field2chk4distr<- field4distrchk ;  fname.prefix <- fname.prefix
#bins=NULL;seedval<-seed
randomisation_ref <- function(q1data.overlap, q2data.overlap, field4mean, field2chk4distr, nPermutations, fname.prefix, seedval, bins=NULL) {
  ##for q1.ref and q2.ref dataset, runs randomisation and returns the mean and median field4mean, saves the obs and rand vals

  #q1data.overlap <- dataframe containing the query1 data that overlaps with the ref
  #i.e ASEs that overlap with the ref database,
  #            avg reads, minFDR
  #q2data.overlap < dataframe being compared against e.g non ASEs that overlap with the ref db
  #e.g.run for ASEs vs nonASEs

  #The q1, q2 dataframes will need to have the following column names
  #"cmp.col" : the "rsid", "chr_pos" fields to be compared

  #field4mean  : : the column name to be used to rank the SNPs/ calculate the mean/median - e.g if spidex "abs_dpsi_max_tissue" or in case of GWAS "neglog10pval" etc
  # field2chk4distr : name of the column to be used to check the distribution of the randomly selected values (in q2) are same as that in the observed ( q1 )
  #               e.g "avgreads" , "MAF"

  #nPermutations : no. of iterations of random selection
  #bins : bins to which loci are to be assigned based on their value in field2chk4distr ,
  #       If not provided and value NULL - code will assign it a binwidth of 2 as long as the field2chk4distr value is  < 200 and remaining in
  #       a single last bin.(used internally for average reads disribution )

  #returns the mean and median of the column in field4mean as a list
  # obs.meanval, ran.meanvals, obs.medianval, ran.medianvals

  #saves the q1data.overlap, ran.from.q2data.overlap , obs.meanval, ran.meanvals, obs.medianval, ran.medianvals
  # in file : <fname.prefix>_obsRandomVals.rda
  #saves the distribution of the mean and median of the obs and random filed4mean
  #in file : <fname.prefix>_DistributionMeanMedian<field4mean>.pdf

  #saves the distribution of the mean nSNPs
  #in file : <fname.prefix>_DistrMeanSNPsPerBin.pdf

  if (!is.null(seedval)){
    set.seed(seedval)
    cat("\n randomisation_ref() Seed set=", seedval, "\n")
  } else {
    cat("\n randomisation_ref()- No seed set \n")
  }
  setDF(q1data.overlap)
  setDF(q2data.overlap)

  #get observed mean  - the q1.ref
  cat("\n mean calculated for : ", field4mean, "\n")
  obs.meanval <- mean( q1data.overlap[ , field4mean] )
  cat("median calculated for: ", field4mean, "\n")
  obs.medianval <- median(  q1data.overlap[ , field4mean])

  #
  setDT(q1data.overlap)
  setDT(q2data.overlap)

  if (is.null( field2chk4distr) ) {
    #no distribution check
    cat("\nRandom SNPs are selected with 'NO' distribution check matching that of the observed.\n")
    nSNPs2Select <- nrow( q1data.overlap )
    perm <- randomisation_ref_noDistributionChk ( q2data.overlap, field4mean, nPermutations, nSNPs2Select )
  } else {
    #if the random are to be selected based on "field2chk4distr" columns distribution of the observed
    cat("\nRandom SNPs with '", field2chk4distr, "' distribution matching that of the observed are selected.\n")

    perm <- randomisation_ref_DistributionChk ( q1data.overlap, q2data.overlap, field4mean, field2chk4distr, nPermutations, bins )
  }
  setDF(q1data.overlap)
  rm(q2data.overlap)

  MEANVAL <- 1  #mean val
  RANDOMHETSOVERLAP <- 2
  MEDIANVAL <- 3  # median val

  #field to compare and field2chk4distr field of the random hets selected
  ran.from.q2data.overlap <-  sapply( perm, function(z) z[RANDOMHETSOVERLAP] )
  ran.from.q2data.overlap <- data.frame( rbindlist(ran.from.q2data.overlap), stringsAsFactors=F )

  #Mean
  ran.meanvals <- unlist (sapply( perm, function(z) z[MEANVAL] ) )
  #Median
  ran.medianvals <-  unlist( sapply( perm, function(z) z[MEDIANVAL] ) )

  #save the perm variable
  q1data.overlap <- q1data.overlap[ , c("cmp.col", field4mean, field2chk4distr) ]
  fname <- paste0(fname.prefix, "_obsRandomVals.rda" )
  save(q1data.overlap, ran.from.q2data.overlap , obs.meanval, ran.meanvals, obs.medianval, ran.medianvals,  file=fname, compress=TRUE)

  cat ("\nOutput dir : ", getwd() , "\nOutput file : ", fname , "\n")

  if (! is.null( field2chk4distr) ) {
    ##plot showing the distribution of (i)mean number of SNPs in randomly selected and observed, (ii) mean no. of unique SNPs
    plotMeanSNPsPerField2Chk4DistrBin ( q1data.overlap, ran.from.q2data.overlap , field2chk4distr, fname.prefix)
  }
  rm(q1data.overlap, ran.from.q2data.overlap )
  gc()

  #Saving the plots for the distribution of the random mean( field4mean) and the obs mean( field4mean) and the median
  plotDistributionMeanMedian( obs.meanval, ran.meanvals, obs.medianval, ran.medianvals, field4mean, fname.prefix )
  return(  list( obs.meanval, ran.meanvals, obs.medianval, ran.medianvals)  )   #return the observed and random means, obs and random median values
}

#function runs the randomisation iteration such that the random are
#selected based on "field2chk4distr" columns distribution of the observed
randomisation_ref_DistributionChk <- function( q1data.overlap, q2data.overlap, field4mean, field2chk4distr, nPermutations, bins ) {
  setDF(q1data.overlap)
  setDF(q2data.overlap)

  nSNPs2Select <- nrow( q1data.overlap )

  distrMatch <- q1data.overlap[ , field2chk4distr]

  #brk.probs <- getProbabilities(distrMatch) # returns a list of breaks and probs
  brk.probs <- getProbabilities(distrMatch, bins ) # returns a list of breaks and probs
  #accessing the 1st element of breaks : brk.probs[[BREAKS]][1]

  #assigning the q2.ref (valid hetSNPs that overlap with ref) an index based on the avg Reads it has using the breaks obtained for the q1.ref (ASEs)
  q2data.overlap <- assignBinIndex (q2data.overlap , brk.probs )   #returns the df with an additonal column containing the bin.id

  nPm <- 1:nPermutations
  cat("\n SNPs2Select = ", nSNPs2Select, "\n")
  cat("\n Randomisation..")
  perm <-  lapply( nPm, function( y, all.hets, nSNPsSel ) {

    #randomly selecting from query2 dataset hets that overlapped to have the same avg read distribution as the q1.ref(ASEs)
    rhets.asedis <- selectRandomValhetWithAseDistrib( all.hets, nSNPsSel, brk.probs)
    #cat("\n", y)

    #get the details
    rhets.asedis <- data.frame( cmp.col = rhets.asedis, stringsAsFactors=F)
    r.het.recs <- merge( all.hets, rhets.asedis, by= "cmp.col", all.y =T )
    rm(rhets.asedis)

    #Get mean and median of field
    val <- r.het.recs[ , field4mean ]
    mean.val <- mean ( val)
    median.val <- median ( val)
    rm(val)
    r.het.recs <- r.het.recs[ , c("cmp.col", field4mean , field2chk4distr, "bin.id") ] #bin id included to generate the plot avg distributions per bin
    #r.het.recs <- r.het.recs[ , c("cmp.col", field4mean , field2chk4distr ) ]
    r.het.recs$perm.no <- y

    vals <- list( mean.val, r.het.recs, median.val )
    rm(mean.val, median.val, r.het.recs)

    if( (y %% 100) == 0 ) {
      gc()
      cat(y, "..")
    }

    return( vals)
  }, all.hets <- q2data.overlap, nSNPsSel <- nSNPs2Select  )
  cat("\n Randomisation end.. \n")
  return(perm)
}

#function runs the randomisation iteration such that the random are randomly selected
#They are NOT selected based on "field2chk4distr" columns distribution of the observed
#q1data.overlap <- q1.ref ; q2data.overlap <- q2.ref; field4mean <- ref.field4mean ; field2chk4distr<- "avg.reads" ;  fname.prefix <- fname.prefix
randomisation_ref_noDistributionChk <- function( q2data.overlap, field4mean, nPermutations, nSNPs2Select ) {

  setDF(q2data.overlap)

  q2SNPs <- q2data.overlap$cmp.col

  nPm <- 1:nPermutations
  cat("\n SNPs2Select = ", nSNPs2Select, "\n")
  cat("\n Randomisation..")
  perm <-  lapply( nPm, function( y, all.hets, nSNPsSel ) {

    #randomly selecting from query2 dataset hets that overlapped
    #rhets.asedis <- sample( all.hets, nSNPsSel, SAMPLEREPEAT )
    rhets.asedis <- sample( all.hets, nSNPsSel, replace=RSID.REPEAT )
    #cat("\n", y)

    #get the details
    rhets.asedis <- data.frame( cmp.col = rhets.asedis, stringsAsFactors=F)
    r.het.recs <- merge( q2data.overlap, rhets.asedis, by= "cmp.col", all.y =T )
    rm(rhets.asedis)

    #Get mean and median of field
    val <- r.het.recs[ , field4mean ]
    mean.val <- mean ( val)
    median.val <- median ( val)
    rm(val)
    #r.het.recs <- r.het.recs[ , c("cmp.col", field4mean , field2chk4distr) ]
    r.het.recs <- r.het.recs[ , c("cmp.col", field4mean ) ] #as no distribution check
    r.het.recs$perm.no <- y
    vals <- list( mean.val, r.het.recs, median.val )
    rm(mean.val, median.val, r.het.recs)
    if( (y %% 100) == 0 ) {
      gc()
      cat(y, "..")
    }

    return( vals)
  }, all.hets <- q2SNPs, nSNPsSel <- nSNPs2Select  )
  cat("\n Randomisation end.. \n")
  return(perm)
}


plotDistributionMeanMedian <- function ( obs.meanval, ran.meanvals, obs.medianval, ran.medianvals, field4mean, fname.prefix )    {
  #generates and saves plots to
  #<fname.prefix>_DistributionMeanMedian<field4mean>.pdf in the current directory

  fname <- paste0(fname.prefix, "_DistributionMeanMedian", field4mean, ".pdf" )
  pdf(file=fname)

  par(mfrow = c(1,1))

  min.val <- min( ran.meanvals,  obs.meanval )
  max.val <- max( ran.meanvals,  obs.meanval )
  x.ax.lim <- c(  min.val - 0.1, max.val +0.1  )
  title.size <- 0.75

  title <- paste0(fname.prefix, "\n:Distribution of mean ", field4mean ,"-Random\n Observed mean ", field4mean ,"=", obs.meanval )
  hist(ran.meanvals, xlim=x.ax.lim , xlab= paste0("mean( ", field4mean ," )"), main= title, cex.main=title.size )
  abline( v=obs.meanval, col="red", lty=2 )
  text(  obs.meanval  , -5, "ASE", 1, col="red", cex=0.6)

  min.val <- min( ran.medianvals,  obs.medianval )
  max.val <- max( ran.medianvals,  obs.medianval )
  x.ax.lim <- c(  min.val - 0.1, max.val +0.1  )
  title <- paste0(fname.prefix, "\n:Distribution of median ", field4mean ,"-Random \n Observed median ", field4mean ,"=", obs.medianval )
  hist(ran.medianvals, xlim=x.ax.lim , xlab= paste0("median( ", field4mean ," )" ), main= title, cex.main=title.size)
  abline( v=obs.medianval, col="red", lty=2 )
  text(  obs.medianval  , -5, "ASE", 1, col="red", cex=0.6)

  dev.off()
  cat ("\nPlots for distribution of the mean & median - obs & random : ", fname , "\n")
  return(T)
  #end saving plot
}



#vals <- distrMatch
getProbabilities <- function( vals, bins) {
  #get the probabilities for the vals average Reads ( depending on field4chkdistr - values for the same will be passed from the calling function),
  #function returns a list of breaks and corresponding probabilities
  #e.g. breaks=0, 50, 100 and probabilities=  0.6815806246, 0.1444231995, 0.0557042702
  #0.6815806246 is the prob of avg reads >0 and <=50 = 0.6815806246

  #bins : if provided (not NULL), this is used for the breaks and to get the probabilities.
  #           if NULL, then BIN.MODE is checked - default is for AVGREADS, can have separate code added if want to set the bins via code for e.g. MAF

  #get the breaks
  if( ! is.null(bins)) {
    brk <- bins
  } else {
      if (BIN.MODE == "AVGREADS"){
          brk <- getBreaksBin( vals)
      }
  }


  h <- hist(vals, breaks=brk, plot=F)
  #breaks <- hist(vals, breaks=brk, plot=F)$breaks
  #counts <- hist(vals, breaks=brk, plot=F)$counts
  #get the probabilities
  prob <- h$counts/sum(h$counts)
  return(list( h$breaks,prob) )
}

getBreaksBin <- function ( vals) {
  #function returns the breaks
  ##setting the bins of 200 width
  #many outliers therefore initial of binwidth as per value in BINWIDTH1 where vals (avgReads) <= 200
  #and remaining in a single bin

  binwidth1 <- BINWIDTH1  #2
  to <- BINWIDTH1_END #commented to make it generic twin data. to <- 200
  #to <- 25 *binwidth1
  brk.set1 <- seq(0, to, by= binwidth1)
  #to plot/check the hist
  #l1 <- vals[ which(vals <=200) ]
  #hist(l1, breaks=brk.set1, main="Hist Distribution of ASEs with avg reads (<=200) \nin PUTM that overlap with the ref dataset")

  binwidth2add <- 1000
  #from <- to + binwidth2
  #brk.set2 <- seq(from, max(vals)+ binwidth2, by= binwidth2)
  brk.set2 <- max(vals) + binwidth2add
  #l1 <- vals[ which(vals >200) ]
  #hist(l1, breaks=brk.set2, main="Hist Distribution of ASEs with avg reads (>200) \nin PUTM that overlap with the ref dataset")
  brk <- c(brk.set1, brk.set2)
  rm(brk.set1, brk.set2, binwidth1, binwidth2add, to)
  return(brk)
}

BREAKS <-1
PROBS <-2

assignBinIndex <- function (df , breaks.probs ) {
  #assign a unique index/binid to the valid het based on the bin its avg reads fall into
  #returns the df with an additonal column containing the bin.id

  #df <- q2data.overlap; breaks.probs <- brk.probs
  lbrk <- length( breaks.probs[[BREAKS]] )-1   #length of the breaks the breaks start with 0 i.e 0,2,4...
  #lbrk <- 6
  cts <- 1:lbrk
  df$bin.id <- NA

  #retval <- lapply( cts, function( i ) {
  for( i in cts){
    #i<-3
    from <- breaks.probs[[BREAKS]][i]
    to <- breaks.probs[[BREAKS]][ i+1]
    #pr <- breaks.probs[[PROBS]][ i]
    r <- which(df$avg.reads > from & df$avg.reads <= to )
    df[ r, "bin.id"] <- i
    rm(from, to, r )
    #       return(i)
  }
  #)
  return(df)
}

selectRandomValhetWithAseDistrib <- function( valid, nSNPs, breaks.probs) {
  #run the sample( vector of the unique bin ids, nSNPs=ASE SNPs that overlap, repeat=T, prob list from the ASEs)
  #sample
  # breaks.probs <- brk.probs; nSNPs<- nSNPs2Select; valid <- all.hets
  ids <- 1: ( length( breaks.probs[[BREAKS]] ) -1 )
  r.ids <- sample(x=ids , size=nSNPs, replace=TRUE, prob=breaks.probs[[PROBS]] )   # bin ids selected based on the probabilities
  r.id.freq <- data.frame( table(r.ids), stringsAsFactors=F)      # freq of the bin ids selected
  rm(r.ids )
  i <- sapply(r.id.freq, is.factor)
  r.id.freq[i] <- lapply(r.id.freq[i], as.character)
  r.id.freq$r.ids <- as.numeric( r.id.freq$r.ids)
  rm(i)


  #Randomly selecting the nvalid hetSNPs ==freq for the bin id using the sample( rsid with the unique id = bin.id, repeat=F, no. of hets to return=freq of bin id ) default distribution.
  nRows <- 1: nrow(r.id.freq)
  ran.hets <- lapply( nRows, function(i ) {
    binid <- r.id.freq[i, "r.ids" ]
    nhets2Select <- r.id.freq[i, "Freq" ]
    rsids.lst <- valid[ which(valid$bin.id == binid), "cmp.col" ]
    #cat(" ", i)
    r.hets <- sample(x=rsids.lst , size=nhets2Select, replace=RSID.REPEAT )
    rm(binid, nhets2Select, rsids.lst )
    return( r.hets)
  } )
  ran.hets <-  unlist( ran.hets) # PUTM-PD = 7845;
  rm(nRows, r.id.freq)
  return( ran.hets)
}


############### Functions - post randomisation function ###############
source(paste0(BASEDIR,"/scripts/randInteRandomisation_Zscorefunctions.R"))




