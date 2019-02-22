#' @import data.table
#' @import ggplot2
#' @import gridExtra
NULL

#Global variables

eraseRsidRepeat <-TRUE      #Variable to set sampling with repeats True or False

#setting for bins
eraseBinwidth <- NULL # 2
eraseBinwidthEnd <- 200

eraseBinMode <- "AVGREADS"
# In the code, if the bin list is provided the use that , else if AVGREADS (set as default)
# then use binwidth 2 (eraseBinwidth) for avgreads  <200 (eraseBinwidthEnd) and remaining in 1 bin
# can add on if internally want to set for MAF instead of passing the bins

#End -setting for bins

eraseBreaks <-1
eraseProbs <-2

#Global variables end


#' Runs the randomization step for Single and multiple annotation
#' enrichment(SAE, MAE respectively).
#' If the mode is SAE, generates the p-value else zscores.
#'
#' The randomization step is run. This is followed by the p-value
#' calulation if the mode is SAE. If the mode is not SAE, assumes
#' MAE and generates the z-scores.
#' The default binwidth is set to 2 (see bins parameter below)
#' The following output files are generated for randomization
#' and saved in the current directory
#'    <outFilePrefix>_obsRandomVals.rda :
#'    <outFilePrefix>_DistributionMean<colname_rankSNPann>.pdf : saves the
#'        distribution of the mean observed and random values
#'    <outFilePrefix>_DistrMeanSNPsPerBin.pdf
#' if mode is SAE
#'    <outFilePrefix>_obsRandomValsZscorePval.rda - contains the zscores for
#'    the observed, randomly selected  values and the pnorm pvalue
#'    <outFilePrefix>_DistributionPvalZscore.pdf
#'  else
#'   <outFilePrefix>_zscoreNmlDistr.rda
#'   <outFilePrefix>_zscoreDistr.pdf

#'    are also generated
#'
#' @param df_sigASE_SNPann  -dataframe containing the significant
#' ASEs that overlap with the SNP annotation, will need to have
#' the column name "cmp.col" (containing e.g. "rsid", "chr_pos"
#' values used to find the interesection with the snp annotation
#' dataset) and those passed via the parameters
#' colname_rankSNPann and colname_chk4distr
#' @param df_nonASE_SNPann -dataframe containing the non significant
#' ASEs that overlap with the SNP annotation, will need to have
#' the column name "cmp.col" and those passed via the parameters
#' colname_rankSNPann and colname_chk4distr
#' @param colname_rankSNPann the column name in the above two
#' data frames with the transformed SNP score,
#' to be used to rank the SNPs/ calculate the mean
#' e.g. in case of GWAS p value can be transformed to
#'  "neglog10pval" containing -log10(p)
#' @param colname_chk4distr name of the column to be used to check the
#' distribution of the randomly selected values (in df_nonASE_SNPann)
#' are same as that in the df_sigASE_SNPann e.g. "averageReads"
#' @param outFilePrefix The names of all the output files
#' generated will be assigned this prefix
#' @param nIterations number of random selection iterations.
#' Default value of 10000
#' @param binwidth Bindwidth to be used. default = 2
#' #' @param mode  - default value of 'SAE', calculates the p-value
#' else for MAE will transform into z-score
#'
#' @param seedValue The seed value to be set. If NULL then no seed
#' is set.
#' @param bins bins to which loci are to be assigned based on their
#' value in colname_chk4distr.
#' If not provided and value NULL - code will assign it a binwidth
#' of 2 ( in variable 'eraseBinwidth') as long as the
#' colname_chk4distr value is  < 200 (in variable 'eraseBinwidthEnd')
#' and the remaining are placed in a single last bin.
#' The defaults can be changed by assigning the
#' variables 'eraseBinwidth' and 'eraseBinwidthEnd' new values.

#' @return the pnorm p-value if mode is 'SAE',
#' transformed zscore rda file name if mode is 'MAE'
#' @export
#'
#' @examples
randomization <- function(df_sigASE_SNPann,
                          df_nonASE_SNPann,
                          colname_rankSNPann,
                          colname_chk4distr,
                          outFilePrefix,
                          nIterations=10000 ,
                          binwidth=2, #eraseBinwidth,
                          mode="SAE",
                          seedValue=NULL,
                          bins=NULL ){

  cat("\n\nERASE mode-", mode, "\nRandomization step start-" )
  date()

  #if( binwidth != eraseBinwidth) {
    eraseBinwidth <<- binwidth
  #}

  obs.ran <- randomisation_ref( df_sigASE_SNPann, df_nonASE_SNPann,
                                colname_rankSNPann, colname_chk4distr,
                                nIterations, outFilePrefix, seedValue , bins)
  date()
  cat("\nRandomization step end \n" )
  rm(obs.ran)
  gc()

  obsranRda <- paste0(outFilePrefix, "_obsRandomVals.rda")
  # if mode= SAE, calculating pval
  #if MAE, generating z-scores
  if(mode=="SAE"){
    cat("\nCalculating p-value :")
    outdir <- getwd()
    ret_val <- getPvalUsingZscore(obsranRda, outdir)
    cat("\n\nPval :", ret_val ,"\n")
    gc()

  } else {
    cat("\nGenerating zscores for: ", obsranRda, "\n")
    outdir <- getwd()
    ret_val <- transform2StdNmlDistr (obsranRda)
  }
  return(ret_val)
  }



# randomisation ref - calls function based on if thhe distribution check is to be done on the randomly selected
#set of SNPs or if no check is to be done
randomisation_ref <- function(q1data.overlap, q2data.overlap, field4mean, field2chk4distr, nPermutations, fname.prefix, seedval=NULL, bins=NULL) {

  #q1data.overlap <- dataframe containing the query1 data that overlaps with the SNP annotation
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
  #in file : <fname.prefix>_DistributionMean<field4mean>.pdf

  #saves the distribution of the mean nSNPs
  #in file : <fname.prefix>_DistrMeanSNPsPerBin.pdf

  if (!is.null(seedval)){
    set.seed(seedval)
    cat("\nSeed set=", seedval, "\n")
  } else {
    cat("\nNo seed set \n")
  }
  cat("Bin width = ", eraseBinwidth)

  setDF(q1data.overlap)
  setDF(q2data.overlap)

  #get observed mean  - the q1.ref
  #cat("\n mean calculated for : ", field4mean, "\n")
  obs.meanval <- mean( q1data.overlap[ , field4mean] )
  #cat("median calculated for: ", field4mean, "\n")
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

  cat ("\nOutput folder : ", getwd() , "\nOutput file : ", fname , "\n")

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
  #accessing the 1st element of breaks : brk.probs[[eraseBreaks]][1]

  #assigning the q2.ref (valid hetSNPs that overlap with ref) an index based on the avg Reads it has using the breaks obtained for the q1.ref (ASEs)
  q2data.overlap <- assignBinIndex (q2data.overlap , brk.probs, field2chk4distr )   #returns the df with an additonal column containing the bin.id

  nPm <- 1:nPermutations
  cat("\nSNPs2Select = ", nSNPs2Select, "\n")
  cat("\nRandomisation..")
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
  cat("\nRandomisation end.. \n")
  return(perm)
}

#function runs the randomisation iteration such that the random are randomly selected
#They are NOT selected based on "field2chk4distr" columns distribution of the observed
randomisation_ref_noDistributionChk <- function( q2data.overlap, field4mean, nPermutations, nSNPs2Select ) {

  setDF(q2data.overlap)

  q2SNPs <- q2data.overlap$cmp.col

  nPm <- 1:nPermutations
  cat("\n SNPs2Select = ", nSNPs2Select, "\n")
  cat("\n Randomisation..")
  perm <-  lapply( nPm, function( y, all.hets, nSNPsSel ) {

    #randomly selecting from query2 dataset hets that overlapped
    #rhets.asedis <- sample( all.hets, nSNPsSel, SAMPLEREPEAT )
    rhets.asedis <- sample( all.hets, nSNPsSel, replace=eraseRsidRepeat )
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



getProbabilities <- function( vals, bins) {
  #get the probabilities for the vals average Reads ( depending on field4chkdistr - values for the same will be passed from the calling function),
  #function returns a list of breaks and corresponding probabilities
  #e.g. breaks=0, 50, 100 and probabilities=  0.6815806246, 0.1444231995, 0.0557042702
  #0.6815806246 is the prob of avg reads >0 and <=50 = 0.6815806246

  #bins : if provided (not NULL), this is used for the breaks and to get the probabilities.
  #           if NULL, then eraseBinMode is checked - default is for AVGREADS, can have separate code added if want to set the bins via code for e.g. MAF

  #get the breaks
  if( ! is.null(bins)) {
    brk <- bins
  } else {
    if (eraseBinMode == "AVGREADS"){
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
  ##setting the bins
  ##initial of binwidth as per value in eraseBinwidth where vals (avgReads) <= 200
  #and remaining in a single bin

  binwidth1 <- eraseBinwidth
  to <- eraseBinwidthEnd
  brk.set1 <- seq(0, to, by= binwidth1)

  binwidth2add <- 1000
  brk.set2 <- max(vals) + binwidth2add
  brk <- c(brk.set1, brk.set2)
  rm(brk.set1, brk.set2, binwidth1, binwidth2add, to)
  return(brk)
}


assignBinIndex <- function (df , breaks.probs, colNameDistrChk ) {
  #assign a unique index/binid to the valid het based on the bin its avg reads fall into
  #returns the df with an additonal column containing the bin.id

  lbrk <- length( breaks.probs[[eraseBreaks]] )-1   #length of the breaks the breaks start with 0 i.e 0,2,4...
  cts <- 1:lbrk
  df$bin.id <- NA
  names(df)[ grep(colNameDistrChk, names(df), fixed = T)] <- "avg.reads"

  for( i in cts){
    from <- breaks.probs[[eraseBreaks]][i]
    to <- breaks.probs[[eraseBreaks]][ i+1]
    r <- which(df$avg.reads > from & df$avg.reads <= to )
    df[ r, "bin.id"] <- i
    rm(from, to, r )
  }
  names(df)[ grep( "avg.reads", names(df), fixed = T)] <- colNameDistrChk

  return(df)
}

selectRandomValhetWithAseDistrib <- function( valid, nSNPs, breaks.probs) {
  #run the sample( vector of the unique bin ids, nSNPs=ASE SNPs that overlap, repeat=T, prob list from the ASEs)
  #sample
  # breaks.probs <- brk.probs; nSNPs<- nSNPs2Select; valid <- all.hets
  ids <- 1: ( length( breaks.probs[[eraseBreaks]] ) -1 )
  r.ids <- sample(x=ids , size=nSNPs, replace=TRUE, prob=breaks.probs[[eraseProbs]] )   # bin ids selected based on the probabilities
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
    if(nhets2Select >0 & length(rsids.lst)==0){
      cat("\n BinWidth = ", eraseBinwidth, "\n")
      stop("Please check bin width. No SNPs in bin to randomly select from.")
    }
    r.hets <- sample(x=rsids.lst , size=nhets2Select, replace=eraseRsidRepeat )
    rm(binid, nhets2Select, rsids.lst )
    return( r.hets)
  } )
  ran.hets <-  unlist( ran.hets) # PUTM-PD = 7845;
  rm(nRows, r.id.freq)
  return( ran.hets)
}


#functions to calculate the pvalue using zscores
#-----------------------------------------------
OBSZ <- 1  #obs zscore
RANZ <- 2 #ran zscore
PNORM <- 3

calculateRanObsZscoresPnorm <- function( obs.val, ran.vals, nOverlaps=NA ) {
  #obs.val = observed mean value
  #ran.vals= list of random means
  #nOverlaps= no. of ASE SNPs that overlapped with the ref/Source [ no. of randomly selected SNPs in each iteration of the bootstrap]
  #           Not in use currently
  #function returns a list of the observed,random zscores, pnorm value

  mean.ranvals <- mean( ran.vals)

  sd.ranvals <- sd( ran.vals)

  #Calculate the z-score for each permutation:
  nPm <- 1: length(ran.vals)
  ran.zscores <-  lapply( nPm, function( y, all.ran ) {
    p <- all.ran[y]
    #Calculate the z-score
    z <- ( p - mean.ranvals)/ (sd.ranvals)
    return( z)
  }, all.ran <- ran.vals  )

  ran.zscores <- unlist (ran.zscores)

    #Calculate the z-score for the observed overlap
  obs.zscore <- ( obs.val - mean.ranvals)/ ( sd.ranvals )

  #calculating the pnorm
  pnorm.val <- pnorm( obs.zscore,lower.tail=F)

  return ( list(obs.zscore, ran.zscores, pnorm.val ))
}


getPval <- function(obs.meanval, ran.meanvals, nOverlaps, fname.prefix ) {
  #nOverlaps = no. of category SNPs that overlap with the source i.e. reference e.g. ASE GWAS overlap
  #fname.prefix = used for the plot title and plot saved as <fname.prefix>_DistributionPvalZscore.pdf
  #               and rda file for the zscores and pval calculated

  #calls the functions to calculate the zscores and pnorm
  # return the pval ( generated by pnorm)
  #Also saves the plot the distribution of the observed and ran zscores

  zscores <- calculateRanObsZscoresPnorm( obs.meanval, ran.meanvals, nOverlaps )
  obs.zscores <-  unlist(zscores[OBSZ] )
  ran.zscores <- unlist ( zscores[RANZ] )
  pnorm.pval <- unlist ( zscores[PNORM] )

  fname <- paste0(fname.prefix, "_DistributionPvalZscore.pdf" )
  title <- paste0(fname.prefix, ":Distribution of Zscores-Random \n Observed Zscore =", obs.zscores )
  plotDistributionVals( obs.zscores, ran.zscores, fname, title )
  cat("\n Zscore distribution plot : ", fname, "\n")

  #save the zscores and pnorm
  fname <- paste0(fname.prefix, "_obsRandomValsZscorePval.rda" )

  save( obs.zscores, ran.zscores, pnorm.pval, file=fname, compress=TRUE)
  cat("\n Zscores, p-value saved in : ", fname, "\n")

  return ( pnorm.pval )

}

getPvalUsingZscore <- function( obsRandomValsrda, outdir) {
  #calls the function to calculate the zscore and the pnorm using the observed and random bootstrap values in the rda file.
  # saves the zscore, pnorm calculated.

  #returns the pnorm pval

  load(obsRandomValsrda)
  fname.prefix <- gsub("_obsRandomVals.rda", "", basename(obsRandomValsrda) )
  currdir <- getwd()
  cat("\n rdafile= ", obsRandomValsrda, "\n")
  cat("Outfut folder= ", outdir , "\n")
  #cat("Output file prefix= ", fname.prefix, "\n")
  setwd(outdir)
  pval <- getPval(obs.meanval, ran.meanvals, nrow(q1data.overlap), fname.prefix )
  setwd(currdir)
  return( pval)
}


#seed for the random number generator from a character string
# if no char string is passed, "RANDOMIZATION" is used as default
getSeed4chrStr<- function(charval="RANDOMIZATION"){
  cat("\nSeed string = ", charval, "\n")
  return( char2seed(charval, set=F) )
}


#Transform -MAE
transform2StdNmlDistr <- function (rda){
  #function loads the randomisation file with the obs and ran mean and median values.
  # calls the function to generate the zscores for the std nml distribution using the mean randomised values.

  #read randomised file to be transformed
  load(rda)   #obs.meanval, ran.meanvals, obs.medianval, ran.medianvals

  fname.base <- gsub( ".rda", "", rda)
  fname.base <- basename(fname.base)

  #getting the std normalised distribution using the medians
  #fname.prefix <- paste0( fname.base, "Median")
  #out.rda.median <- zscoreStdNmlDistr (obs.medianval, ran.medianvals, fname.prefix)

  #getting the std normalised distribution using the means
  fname.prefix <- paste0( fname.base, "Mean")
  out.rda <- zscoreStdNmlDistr (obs.meanval, ran.meanvals, fname.prefix)

  cat ("\n Zscore generation complete \n")
  return(out.rda)

}

#calculates the zscores to transform the distribution to a standard normal distribution
# also plots the ran and obs zscore distribution
#obs<- obs.medianval ; ran <- ran.medianvals;  file.prefix <-  fname.prefix
zscoreStdNmlDistr  <- function(obs, ran, file.prefix){
  #function calculates the Z-score for the observed and the random set of values
  #obs : obs value obtained from the randomisation step,
  #ran : random values generated from the randomisation step
  # saves the obs and random zscores to a rda file  <file.prefix>_zscoreNmlDistr.rda
  # returns the rda file name
  #distribution of the zscores saved to : <file.prefix>_zscoreDistr.pdf

  #z-score calculation:
  # v1, v2, v3.. vn : n=no. of randomisations
  #Z-score = ( val -  X )/(sd)
  #val e.g. v1
  #X = mean(v1, v2, v3.. vn)
  #sd =  sd( v1, v2, v3.. vn)

  #calculate the zscores for the random values
  nPm <- 1:length( ran)
  mu <- mean( ran )
  sd <- sd( ran )
  ran.zscores <-  lapply( nPm, function( y, all.vals ) {
    val <- all.vals[y]    # val for this randomisation
    #Calculate the z-score
    z <- ( val - mu)/ ( sd )
    return( z)
  }, all.vals <- ran  )
  ran.zscores <- unlist( ran.zscores)

  #Calculate the z-score for the observed value
  obs.zscore <- (obs - mu)/sd

  fname <- paste0(file.prefix, "_zscoreNmlDistr.rda" )

  save(obs.zscore, ran.zscores, file=fname, compress=TRUE)
  cat ("\nOutput dir : ", getwd() , "\nOutput file : ", fname , "\n")

  p.fname <- paste0(file.prefix, "_zscoreDistr.pdf" )
  title <- gsub(".pdf" , "", p.fname)
  title <- paste0( title, "\n: observed value = ", obs.zscore )
  plotDistributionVals (obs.zscore, ran.zscores, p.fname, title)

  cat ("\nZscore distribution plot saved to : ", p.fname , "\n")

  return ( fname )
}


## Functions for the plots
#--------------------------



#Functions


#plot showing the distribution of mean number of SNPs per field2check bin e.g. avg.reads
#of the randomly selected SNPs [bins correspond to avg read depth]
#across the 10K randomisations and obs

#Plot 2 - mean proportion of unique SNPs per bin
#returns the name of the file generated

plotMeanSNPsPerField2Chk4DistrBin <- function( obs,ran , field2chk4distr, fname.prefix, bins=NULL ) {

  names(obs)[ grep(field2chk4distr, names(obs), fixed = T)] <- "avg.reads"
  names(ran)[ grep(field2chk4distr, names(ran), fixed = T)] <- "avg.reads"


  setDT(obs); setDT(ran)
  #setnames(obs, old = field2chk4distr, new = "avg.reads")

  #setnames(ran, old = field2chk4distr, new = "avg.reads")


  #get the avg no. of SNPs per bin across perm
  totPerPermBin <- ran[ , .(nSNPs.perbin.perm=length(cmp.col)), by="perm.no,bin.id"]
  avgPerBin <- totPerPermBin[ , .(nSNPs=mean(nSNPs.perbin.perm)), by="bin.id"]
  avgPerBin$type <- "Mean.nSNPs.perBin.acrossRand"
  setorder(avgPerBin, bin.id)

  #get the mean avg Reads per bin
  meanAvgReadsPerPermBin <- ran[ , .(meanAvgReads.perbin.perm=mean(avg.reads)), by="perm.no,bin.id"]
  meanAvgReadsPerBin <- meanAvgReadsPerPermBin[ , .(meanAvgReads=mean(meanAvgReads.perbin.perm)), by="bin.id"]
  meanAvgReadsPerBin$type <- "Mean.AvgReads.perBin.acrossRand"
  setorder(meanAvgReadsPerBin, bin.id)

  #end get mean avg reads

  #mean proportion of unique SNPs per bin
  uniqPerPermBin <- ran[ , .(nUniqSNPs.perbin.perm=length(unique(cmp.col))), by="perm.no,bin.id"]

  setkeyv(uniqPerPermBin, c("perm.no", "bin.id"))
  setkeyv(totPerPermBin, c("perm.no", "bin.id"))
  res <- totPerPermBin[uniqPerPermBin ]
  res$prop.uniq.SNPs <- res$nUniqSNPs.perbin.perm /res$nSNPs.perbin.perm
  avgPropUniqPerBin <- res[ , .(prop.uniq.SNPs.mean=mean(prop.uniq.SNPs)), by="bin.id"]

  rm(res)

  #obs
  val<- obs[ ,.SD, .SDcol="avg.reads"]
  if( ! is.null(bins)) {
    brk <- bins
  } else {
    if (eraseBinMode == "AVGREADS"){
      brk <- getBreaksBin(val )
    }
  }

  h <- hist(val$avg.reads, breaks=brk, plot=F)
  obsPerBin <- h$counts
  obsPerBin <- data.table( nSNPs = obsPerBin)
  obsPerBin$bin.id <- as.numeric(rownames(obsPerBin))
  setcolorder(obsPerBin, c( "bin.id", "nSNPs" ))
  obsPerBin$type <- "obs.nSNPs.perBin"


  setDF(val)
  #passing the bin id for the avg reads in the observed
  brk.probs <- getProbabilities(val$avg.reads, bins )
  obsmeanAvgReadsPerBin <- assignBinIndex (val , brk.probs, field2chk4distr ) #assign the bin ids for the avg reads using the breaks vals
  #get the avg per bin
  names(obsmeanAvgReadsPerBin)[ grep(field2chk4distr, names(obsmeanAvgReadsPerBin), fixed = T)] <- "avg.reads"

  setDT(obsmeanAvgReadsPerBin)
  obsmeanAvgReadsPerBin <- obsmeanAvgReadsPerBin[ , .(meanAvgReads=mean(avg.reads)), by="bin.id"]
  setorder(obsmeanAvgReadsPerBin, bin.id)
  obsmeanAvgReadsPerBin$type <- "obs.avgReads.perBin"
  #Kolmogorov Smirnov test - avgreads
  pval_avgreads <- ks.test( obsmeanAvgReadsPerBin$meanAvgReads, meanAvgReadsPerBin$meanAvgReads)$p.value


  obsPerBin <- obsPerBin[ bin.id %in% obsmeanAvgReadsPerBin$bin.id, ] #added as the obserBin has the initial bins with 0 SNPs and affects the ks.test
  pval <- ks.test( obsPerBin$nSNPs, avgPerBin$nSNPs)$p.value


  comb <- rbindlist( list(avgPerBin, obsPerBin))


  #Variables for the plots
  MULTI <- 1

  TITLE_SIZE <- 10*MULTI
  LEGEND_TITLE_SIZE <- 10*MULTI
  LEGEND_TEXT_SIZE <- 8*MULTI
  AXIS_TITLE_SIZE <- 10*MULTI
  AXIS_TEXT_SIZE <- 8*MULTI
  POSITION_DODGE_WIDTH <- 0.5

  #theme_clearBg <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #                       panel.background = element_blank(), axis.line = element_line(colour = "black"))

  theme_legend <- theme(legend.title = element_text(colour="black", size=LEGEND_TITLE_SIZE, face="bold"),
                        legend.text = element_text(colour="black", size=LEGEND_TEXT_SIZE))

  theme_title <- theme(plot.title = element_text(size=TITLE_SIZE, face="bold", hjust = 0.5)) #hjust =0.5 to centre
  theme_subtitle <- theme(plot.subtitle=element_text(size=TITLE_SIZE-1, hjust=0.5, face="italic"))
  #X and Y axes labels, text
  theme_axes <- theme(axis.title = element_text(face="bold", size=AXIS_TITLE_SIZE),
                      axis.text  = element_text( vjust=0.5, size=AXIS_TEXT_SIZE) )

  fill_legend <- scale_fill_manual(values=c("turquoise4", "yellow3"))
  fill_legend_meanavgRds <- scale_fill_manual(values=c("seagreen4", "slateblue3"))

  #End variables

  stitle <- basename(fname.prefix)
  title <- paste0( "Distribution of mean #SNPs per bin (", field2chk4distr, ") of \n the randomly selected and the observed")
  title <- paste0(title, "\nKolmogorov Smirnov test p-value=", pval )

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

  ggsave( filename = outfile, marrangeGrob(grobs = plots, nrow=1, ncol=1), width=11, height=8.5)

  return(outfile)
}


plotDistributionMeanMedian <- function ( obs.meanval, ran.meanvals, obs.medianval, ran.medianvals, field4mean, fname.prefix )    {
  #generates and saves plots to
  #<fname.prefix>_DistributionMean<field4mean>.pdf in the current directory

  fname <- paste0(fname.prefix, "_DistributionMean", field4mean, ".pdf" )
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

  # min.val <- min( ran.medianvals,  obs.medianval )
  # max.val <- max( ran.medianvals,  obs.medianval )
  # x.ax.lim <- c(  min.val - 0.1, max.val +0.1  )
  # title <- paste0(fname.prefix, "\n:Distribution of median ", field4mean ,"-Random \n Observed median ", field4mean ,"=", obs.medianval )
  # hist(ran.medianvals, xlim=x.ax.lim , xlab= paste0("median( ", field4mean ," )" ), main= title, cex.main=title.size)
  # abline( v=obs.medianval, col="red", lty=2 )
  # text(  obs.medianval  , -5, "ASE", 1, col="red", cex=0.6)

  dev.off()
  cat ("\nPlots for distribution of the mean - observed & random values: ", fname , "\n")
  return(T)
  #end saving plot
}

plotDistributionVals <- function( obs.val, ran.vals, fname, plot.title ) {
  #function to plot the distribution of the observed and random values provided
  # fname = pdf file name to save the plot
  #plot.title = title to be assigned to the plot

  min.val <- min( ran.vals,  obs.val )
  max.val <- max( ran.vals,  obs.val )
  x.ax.lim <- c(  min.val - 0.1, max.val +0.1  )
  title.size <- 0.75

  pdf(file=fname)

  hist(ran.vals, xlim=x.ax.lim , xlab= "Randomly selected values", main= plot.title, cex.main=title.size  )
  abline( v= obs.val, col="red", lty=2 )
  text( obs.val  , -5, "Observed\nvalue", 1, col="red", cex=0.6)

  dev.off()

  return(T)
}

