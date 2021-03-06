% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intersection.R
\name{getIntersection}
\alias{getIntersection}
\title{Get the SNP set common to/intersection of the ASE calling and
SNP annotation datasets}
\usage{
getIntersection(df_ase, df_snpAnn, snpAnnType = "GWAS",
  snpAnnCol2Sel = "p")
}
\arguments{
\item{df_ase}{- a data frame object containing the SNPs from the ASE
calling dataset with a column named 'cmp.col' e.g. containing
the rsid, chr:pos values to be compared with the SNP Annotation data
It should also have the following columns:
1. with values to be accounted for
during randomization e.g. the average read depth and
2. with values to assess the ASE significance e.g. FDR, p-value}

\item{df_snpAnn}{-  data frame object with the SNP annotation data.
It must have the column name 'cmp.col' - see df_ase above and a column
to rank the SNPs by during randomization}

\item{snpAnnType}{- is the code for the SNP annotation default is GWAS
if the snpAnnType is GWAS, then first, the pvals with value 0  are
updated to next.min*1e-100 ;
SNPs with > 1 pval are assigned the min pval}

\item{snpAnnCol2Sel}{- the column name in the SNP annotaion data
to be used to rank the SNPs during randomization- e.g. 'p' for GWAS
default is 'p' i.e. p-value for GWAS , the default snpAnnType}
}
\value{
a dataframe object with SNPs common to the ASE calling and
SNP annotation dataset with the snpAnnCol2Sel added to the columns
present in df_ASE
}
\description{
Get the SNP set common to/intersection of the ASE calling and
SNP annotation datasets
}
