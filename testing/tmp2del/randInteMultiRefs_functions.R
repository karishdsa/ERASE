
library(data.table)


##function to plot Avg read distribution of q1.ref vs q2.ref
library(ggplot2)



plotHistpval <- function(df, fname) {
#df <- toplot
pdf(file=fname)
b.title <- paste0( gsub(".pdf", "", fname) )

ggplot(df, aes(x=neglogpval, fill=hets)) +
  ggtitle(b.title) +
  geom_histogram()+
  facet_grid(hets~.)

dev.off()

}

plotHistAvgReads <- function(df, fname) {
	      #df <- toplot
	 pdf(file=fname)
	b.title <- paste0( gsub(".pdf", "", fname) )
	ggplot(df, aes(x=avg.reads, fill=hets)) +
		ggtitle(b.title) +
		 geom_histogram(binwidth=10)+
		 facet_grid(hets~.)

	title <- paste0(b.title, "\n Zooming ,xlim=200")
	ggplot(df, aes(x=avg.reads, fill=hets)) + xlim(0,200)+
		ggtitle(title) +
		geom_histogram(binwidth=20)+
		facet_grid(hets~.)

	title <- paste0(b.title, "\n ,ASEs")
	df <- df[ which(df$hets=='ASEs'), ]
	ggplot(df, aes(x=avg.reads, fill=hets)) +
		ggtitle(title) +
		geom_histogram(binwidth=20)

	title <- paste0(b.title, "\n ,ASEs, Zooming")
	df <- df[ which(df$hets=='ASEs'), ]
	ggplot(df, aes(x=avg.reads, fill=hets)) + xlim(0,200)+
		ggtitle(title) +
		geom_histogram(binwidth=20)

	dev.off()
	return(T)
}




#file <- gwas.fname
readGWAS<- function(file, file.delim="tab"){
  #functioned reads the gwas file, and returns the df with the pval field to p & rsid field to MarkerName for consistetncy

  if (file.delim =="space") {
    sepa <- " "
  } else {
    sepa <- "\t"
  }

  if( length( grep(".gz$|.zip$", file)) != 0 ) {
    cmd <- paste0("zcat ", file)
   # df <- read.delim( pipe(cmd), as.is=T, sep=sepa, stringsAsFactors=F)   #
    df <- fread( cmd, sep=sepa )
  } else {
   # df <- read.delim(file, as.is=T, sep=sepa, stringsAsFactors=F)   #

    df <- fread( file, sep=sepa )
  }
  setDF(df)
  names(df)[ grep("^PValue$|SCAN.P|Pval|P-val", names(df), ignore.case=T) ] <- "p"
  names(df)[ grep("^Marker$|SNP", names(df)) ] <- "MarkerName"
  return(df)
}



updatepvaleq0tofixedmin <- function(df) {
  #Assign new min if pval==0 to the smallest number handled in R. If it is greater than the min pval after 0 then flag up
  if( length(which(df$p==0) ) >0 ) {
    #cat("in")
    next.min <- min( df[ which(df$p != 0),"p" ] )
    #new.min <- next.min*1e-100

    new.min <- 5e-324 # smallest in R to 0
    cat("\n Next min p-val =", next.min, "\n New min =", new.min , "\n")
    if (next.min < new.min) {
      cat("\n !!!! POSSIBLE ERROR - THe new min p-val set where p=0 is not the smallest")
    }

    df[ which(df$p==0), "p"] <- new.min
  }
  return(df)
}


getHetOverlapRef <- function( df, p.df, field2sel) {
  #the function returns a dataframe with the overlapping SNPs i.e the SNPs in df present in p.df,
  #and the field2sel col e.g pval from the pval.df using the cmp.col
  # In the df returned the fields will be same as df with an additional column of pval (p)

  #p.df <- p.df[ , c("cmp.col", "p") ]
  p.df <- p.df[ , c("cmp.col", field2sel) ]
  overlap  <- merge( df , p.df, by="cmp.col")   #
  df <- overlap

  ##with( overlap,( tapply( p, cmp.col, function(z) { length(z) })) ) # to check nPvals present per SNP from the overlap

  return(df)
}



getIntersectionValRefs <- function(ref1.rda, ref2.rda, out.rda ){
    # for validGWAS overlap and the validSPIDEX /validref2 overlap - returns the valid.GWAS.ref2
  #ref1.rda <- "/SAN/neuroscience/ukbecvc/ase57/results/randomisationInte/overlapFiles/ACROSS_TISS_validhets_SCZ_Overlap.rda" ;
  #ref2.rda <- "/SAN/neuroscience/ukbecvc/ase57/results/randomisationInte/overlapFiles/ACROSS_TISS_validhets_SPIDEX_Overlap.rda" ;

  load(ref1.rda)  #valid.gwas
  valid.gwas <- valid.gwas [ , c("cmp.col", "p", "ref")]
  #names(valid.gwas)[grep( "ref", names(valid.gwas) )]  <- "ref1"

  load(ref2.rda)  #valid.ref - spi
  #names(valid.ref)[grep( "ref", names(valid.gwas) )]  <- "ref1"

  valid.ref <- merge( valid.gwas, valid.ref, by="cmp.col" )
  valid.ref$ref <- paste0(valid.ref$ref.x, "_", valid.ref$ref.y)
  valid.ref$ref.x <- NULL
  valid.ref$ref.y <- NULL

  cat("valid intersection ", unique(valid.ref$ref ), " :\n ASEs = ", nrow(valid.ref[ which(valid.ref$min.FDR < 0.05), ] ))
  cat("\n Non ASEs = ", nrow(valid.ref[ which(valid.ref$min.FDR >= 0.05), ] ))

  save(valid.ref, file=out.rda, compress=TRUE)
  cat( "output dir : ", getwd(), " , File : ", out.rda)
  return(T)
}



#' Returns the records for a GWAS file
#'
#'The function reads the gwas summary file, updates the pvals with value 0 to next.min*1e-100
#'
#' @param fname
#' @param delim
#'
#' @return
#' @export
#'
#' @examples
getGWASrecs <- function(fname, delim="tab" ){

  gwas <- readGWAS(fname,  delim)
  gwas <- updatepvaleq0 (gwas)
  return(gwas)
}


#####################Functions - randomisation MAE ################################


#Transform
#calculates the zscores to transform the ditribution to a standard normal distribution
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
  plotDistributionObsRan (obs.zscore, ran.zscores, p.fname)

  cat ("\nZscore distribution plot saved to : ", p.fname , "\n")

  return ( fname )
}

transform2StdNmlDistr <- function (rda){
  #function loads the randomisation file with the obs and ran mean and median values.
  # calls the function to generate the zscores for the std nml distribution using the mean and median randomised values.

  #read randomised file to be transformed
  load(rda)   #obs.meanval, ran.meanvals, obs.medianval, ran.medianvals



  fname.base <- gsub( ".rda", "", rda)
  fname.base <- basename(fname.base)

  #getting the std normalised distribution using the medians
  fname.prefix <- paste0( fname.base, "Median")
  out.rda.median <- zscoreStdNmlDistr (obs.medianval, ran.medianvals, fname.prefix)

  #getting the std normalised distribution using the means
  fname.prefix <- paste0( fname.base, "Mean")
  out.rda <- zscoreStdNmlDistr (obs.meanval, ran.meanvals, fname.prefix)

  cat ("\n Complete \n")
  return(T)

}



#Plot and saving the distribution of the obs and ran values
# obs <- obs.zscore;  ran <- ran.zscores
plotDistributionObsRan <- function(obs, ran, outfile.pdf){
    #function saves the distribution of the obs and ran values to
    # outfile.pdf : a filename with .pdf extension
    # in the current directory

    pdf(file=outfile.pdf)

    par(mfrow = c(1,1))

    min.val <- min( ran,  obs )
    max.val <- max( ran,  obs )
    x.ax.lim <- c(  min.val - 0.1, max.val +0.1  )
    title.size <- 0.75

    title <- gsub(".pdf" , "", outfile.pdf)
    title <- paste0( title, "\n: obs value = ", obs )
    hist(ran, xlim=x.ax.lim , xlab= "Random vals", main= title, cex.main=title.size )
    abline( v=obs , col="red", lty=2 )
    text(  obs  , -5, "Obs", 1, col="red", cex=0.6)

    dev.off()
#    cat ("\n Distribution of the - obs & random saved to : ", outfile.pdf , "\n")
    return (T)
}



#Integrate
#obs.ref1<-obs.r1.mean; ran.ref1 <-ran.r1.mean ;obs.ref2<- obs.r2.mean; ran.ref2 <- ran.r2.mean; fname.prefix <-fname
integrate <- function(obs.ref1, ran.ref1, obs.ref2, ran.ref2, alpha, nIterations, fname.prefix) {
  # For a  alpha, integrate obs values (obs.ref1, obs.ref2) for the 2 ref datasets and random values for the same - save the integrated vals

  #alpha = measure of confidence in the ref1
  #(1-alpha) = measure of confidence in the ref2

  #integrated.obs = (obs.ref1 * alpha )+ (obs.ref2 (1-alpha) )

  #nIterations = no. of iterations for the pairwise integration of the random values ( ran.ref1, ran.ref2)
  #For the nIter, randomly select an element from ran.ref1 and ran.ref2,
  #integrated.ran =   (ran.ref1 * alpha )+ (ran.ref2 (1-alpha) )

  #integrated values and alpha saved
  #in file : <fname.prefix>_integratedVals.rda

  # return the rda file name the integrated values are saved - as a list of the integrated.obs, integrated.random

  seedstr <- paste0(SEEDSTRING, as.character(alpha*1000))
  seedval <- getSeed4chrStr(seedstr)
  set.seed(seedval)
  cat("\nSeed set :", seedval, "\n")

  int.obs <- (obs.ref1 * alpha )+ (obs.ref2  *(1-alpha) )

  nIter <- 1:nIterations
  nElts2Sel <- 1  # nElements selected is 1 for each iteration
  int.ran <-  lapply( nIter, function( y ) {
    r.r1 <-  sample(x=ran.ref1 , size=nElts2Sel, replace=RANDPAIRINT.REPEAT )
    r.r2 <-  sample(x=ran.ref2 , size=nElts2Sel, replace=RANDPAIRINT.REPEAT )
    i.ran =   (r.r1 * alpha )+ (r.r2 * (1-alpha) )
    return(i.ran)
  } )
  int.ran <-  unlist( int.ran)

  fname <- paste0(fname.prefix, "_integratedVals.rda" )
  save(int.obs, int.ran , alpha,  file=fname, compress=TRUE)

  cat ("\nOutput dir : ", getwd() , "\nOutput file : ", fname , "\n")

  #Saving the distribution of the integrated obs and ran
  p.fname <- paste0(fname.prefix, "_intValDistr.pdf" )
  plotDistributionObsRan (int.obs, int.ran, p.fname)

  cat ("\nIntegrated values distribution plot saved to : ", p.fname , "\n")

  return( fname )
}

runIntegration4alpha <- function ( ref1.rda,ref2.rda, alpha, nIterations, ref1, ref2, tiss) {
  #Run integrate() to integrate the randomisation done with ref1 and ref2
  # for an alpha value

  # ref1 <- "SCZ2018";  ref2 <- "SPIDEX"; tiss <- "PUTM"
  cat("\n Tissue = ", tiss )
  cat( "\n\nIntegrating ref1 :", ref1 , "and ref2 :" , ref2 , "\n")
  cat("\n ref1 = ", ref1.rda , "\n ref2 = ", ref2.rda )
  cat("\n nIterations = ", nIterations, "\n alpha = ", alpha , "\n")

  #<tiss>_ASEsNonASEs_<ref1><ref2>_Randomisation_nIter<nIterations>_alpha<alpha>
  fname.prefix <- paste0(tiss, "_ASEsNonASEs_", ref1, "_", ref2, "_Randomisation",  "_nIter", nIterations, "_alpha", alpha )

  # Ref1 randomisation values
  load(ref1.rda)
  obs.r1 <- obs.zscore
  ran.r1 <- ran.zscores

  rm(obs.zscore, ran.zscores)

  # Ref2 randomisation values
  load(ref2.rda)
  obs.r2 <- obs.zscore
  ran.r2 <- ran.zscores

  rm(obs.zscore, ran.zscores)

  #Running integration
  fname <- fname.prefix
  cat ("\n Integration start : ", date(), "\n")
  retval <- integrate(obs.r1, ran.r1, obs.r2, ran.r2, alpha, nIterations, fname)
  cat ("\n Integration end :",  date(), "\n")

  rm( obs.r1, ran.r1, obs.r2, ran.r2, alpha, nIterations, fname)
  return(retval)
  ##End running integration
}


#pval

calculatePval4Alpharda <- function( rdalst) {
  #function returns a df of pval calculated usign the integrated rda file
  #gets the alpha value stored in the rda

  pvallist <- lapply(rdalst, function(rdafile ) {
    cat(rdafile, "\n")
    load(rdafile)
    vals <- calculateRanObsZscoresPnorm ( int.obs, int.ran, nOverlaps ) # returns the (obs.zscore, ran.zscores, pnorm.val )
    pnorm.pval <- unlist ( vals[PNORM] )
    return( list( alpha, pnorm.pval ))

  })
  pvals <- as.data.frame.matrix( rbindlist(pvallist), stringsAsFactors=F )
  names(pvals) <- c("alpha", "pval")
  return(pvals)
}


calculatePval4TransScores <- function( rdalst) {
  #function returns a df of pval calculated usign the transformed zscores

  pvallist <- lapply(randFiles, function(rdafile ) {
    #rdafile<- randFiles[1]
    cat(rdafile, "\n")
    load(rdafile)
    vals <- calculateRanObsZscoresPnorm ( obs.zscore, ran.zscores ) # returns the (obs.zscore, ran.zscores, pnorm.val )
    pnorm.pval <- unlist ( vals[PNORM] )
    return( list( rdafile, pnorm.pval ))

  })
  pvals <- as.data.frame.matrix( rbindlist(pvallist), stringsAsFactors=F )
  names(pvals) <- c("rda", "pval")
  return(pvals)
}
#####################Functions - Randomisation ################################
#Randomisation function
source(paste0(BASEDIR,"scripts/randInteRandomisation_functions.R"))

