

#' Get the SNP set common to/intersection of the ASE calling and
#' SNP annotation datasets
#'
#' @param df_ase - a data frame object containing the SNPs from the ASE
#' calling dataset with a column named 'cmp.col' e.g. containing
#' the rsid, chr:pos values to be compared with the SNP Annotation data
#' It should also have the following columns:
#' 1. with values to be accounted for
#' during randomization e.g. the average read depth and
#' 2. with values to assess the ASE significance e.g. FDR, p-value
#' @param df_snpAnn -  data frame object with the SNP annotation data.
#' It must have the column name 'cmp.col' - see df_ase above and a column
#' to rank the SNPs by during randomization
#' @param snpAnnType - is the code for the SNP annotation default is GWAS
#'  if the snpAnnType is GWAS, then first, the pvals with value 0  are
#'  updated to next.min*1e-100 ;
#'  SNPs with > 1 pval are assigned the min pval
#' @param snpAnnCol2Sel - the column name in the SNP annotaion data
#' to be used to rank the SNPs during randomization- e.g. 'p' for GWAS
#' default is 'p' i.e. p-value for GWAS , the default snpAnnType
#' @return a dataframe object with SNPs common to the ASE calling and
#' SNP annotation dataset with the snpAnnCol2Sel added to the columns
#' present in df_ASE
#' @export
#'
#' @examples
getIntersection <- function(df_ase, df_snpAnn, snpAnnType="GWAS", snpAnnCol2Sel="p"){
  if (snpAnnType == "GWAS") {
    df_snpAnn <- updatepvaleq0 (df_snpAnn)
    df_snpAnn <- getMinPval ( df_snpAnn)
  }
  df_inter <- getHetOverlapRef(df_ase, df_snpAnn, snpAnnCol2Sel )
  cat("\n No. of SNPs common =", nrow(df_inter ), "\n")

  return(df_inter)
}



#' Get the intersection for MAE -of the ASE calling, SNP annotation1
#'  and SNP annotation2
#'
#'The function takes as input 2 data frames, one with the
#'intersection of the ASE calling dataset and the SNP annotation1
#'and the second the intersection of the ASE calling and SNP
#'annotation2 dataset. These can be obtained by calling the
#'getIntersection(). Each of the data frame objects must have
#'a column named 'cmp.col', e.g. containing
#'the rsid, chr:pos values, that is used to merge the data.
#'
#' @param df_aseSnpAnn1 dataframe object containing the intersection
#' of the ASE calling dataset and SNP annotation1 -the output of
#' getIntersection()
#' @param df_aseSnpAnn2 dataframe object containing the intersection
#' of the ASE calling dataset and SNP annotation2 -the output of
#' getIntersection()
#' @return a dataframe object with SNPs common to df_aseSnpAnn1 and
#' df_aseSnpAnn2
#' @export
#'
#' @examples
getIntersectionMae <- function(df_aseSnpAnn1, df_aseSnpAnn2){

  return( merge( df_aseSnpAnn1, df_aseSnpAnn2, by="cmp.col" ) )

}


updatepvaleq0 <- function(df) {
  #if pval==0, assign the pval as (next min pval)*1e-100
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

