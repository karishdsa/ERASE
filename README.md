# ERASE

ERASE assesses ASE data for enrichment of SNPs in a single or multiple annotation datasets,
with the aim of gaining insight into its specific properties. It does this accounting for 
read depth an important bias in ASE detection. 

ERASE is described in our paper, ERASE: Extended Randomization for assessment of 
multiple annotation enrichment in ASE datasets <link>


## Install
devtools::install_github("kmdsa/ERASE")

## Usage
This section describes the steps to run ERASE examining enrichment in single (SAE) 
and multiple (MAE) annotation datasets.


### ERASE - single annotation enrichment (SAE)
To test if the ASEs are enriched for SNPs in an annotation dataset e.g. a GWAS summary 
dataset.

**1. Intersection :** Get the overlap of the SNPs examined for ASE(ASE dataset) and the SNP annotation dataset using *getIntersection()*. 

   #### Input files  
   * ASE dataset : This file should have all the SNPs that were examined for ASE and the following :  
     - column name 'cmp.col' that contains e.g. the rsid or chr:pos values to compare with that in the SNP annotation dataset
     - column with values to be accounted for during randomization e.g. the average read depth and 
     - column with values to assess the ASE significance e.g. FDR, p-value
   * SNP annotation dataset : e.g. GWAS summary dataset. It should contain the following :
     - column name 'cmp.col' that contains e.g. the rsid or chr:pos values to compare with that in the ASE dataset.
     - column with values ranking the SNPs e.g. p-values

  
   #### Example
   *ase_annotation <- getIntersection(<df_ase>, <df_snpAnn>)*  
   where *<df_ase>* and *<df_snpAnn>* are data frame objects corresponding to the the ASE dataset.  
   The default SNP annotation dataset is 'GWAS' with 'p' as the column header for values ranking its SNPs.  
   
   See help() for more details, optional parameters and their defaults.  
        
**2. Randomization and p-value calculation :** Run the Randomization process with *randomization()*, 
using the overlapping SNPs obtained from Step 1 above, to get the p-value for enrichment. 
  
   #### Example
   *randomization( <df_sigASE_ann>, <df_nonASE_ann>, <colname_rankSNPann>, <colname_chk4distr>, <outFile_prefix>)*    
   where     
   *<df_sigASE_ann>* and <df_nonASE_ann> are data frame objects containing the significant and non significant ASE signals - the output of *getIntersection()*  
   *<colname_rankSNPann>* is the name of the column with the transformed SNP score e.g. in case of GWAS the p-value can be transformed to "neglog10pval" containing -log10(p), that is used to rank the SNPs  
   *<colname_chk4distr>* name of the column containing values to be accounted for during randomization e.g. the average read depth  
   *<outFile_prefix>* the name to be prefixed to all the output files generated.  
   Randomization is run for a default 10000 iterations with bin size of 2. 
  
   See help() for more details, optional parameters, their defaults and output files generated.  



### ERASE - multiple annotation enrichment (MAE)
To test if the ASEs are enriched for SNPs in two annotation datasets e.g. GWAS summary and transcript 
specificity datasets.

**1. Intersection :** Get the SNPs common to the ASE and two annotation datasets by running *getIntersectionMae()*.
The function takes as input the overlap of the SNPs examined for ASE and each of the SNP annotation dataset
separately. This can be obtained by running *getIntersection()* for each SNP annotation and ASE dataset combination.

   #### Example
   *ase_annotation1 <- getIntersection(<df_ase>, <df_snpAnn1>)*  
   *ase_annotation2 <- getIntersection(<df_ase>, <df_snpAnn2>)*  
  
   *ase_ann1_ann2 <- getIntersectionMae(ase_annotation1, ase_annotation2)*  
   
   See help() for more details, optional parameters and their defaults.
  
**2. Randomization and Transformation :** Run the Randomization process with *randomization()* twice, using 
the overlapping SNPs obtained from the previous step, to get the scores for each of the SNP annotations.
*randomization()* called with the parameter 'mode' set to 'MAE' generates the transformed z-scores.
  
   #### Example
   *randomization()* is run using *ase_ann1_ann2* obtained from *getIntersectionMae()*, once using the transformed SNP scores of annotation1 to rank the SNPs as in (a) below  
   
   (a)  *randomization( <sigASE_ann>, <nonASE_ann> , <colname_rankSNPann1>, <colname_chk4distr>, <outFilePrefix_ann1>, eraseMode="MAE" )*  
   
   and the second time using the transformed SNP scores of annotation2 as in (b)  
   
   (b)  *randomization( <sigASE_ann>, <nonASE_ann> , <colname_rankSNPann2>, <colname_chk4distr>, <outFilePrefix_ann2>, eraseMode="MAE" )*  
   
   where  
   *<sigASE_ann>*, *<nonASE_ann>* are data frame objects containing the significant and non significant ASE signals - the output of *getIntersectionMae()*  
   *<colname_rankSNPann1>*,  *<colname_rankSNPann2>*  are the names of the columns with the transformed SNP scores for annotation 1 and annotation2 respectively that is used to rank the SNPs.  
   *<colname_chk4distr>* name of the column containing values to be accounted for during randomization e.g. the average read depth  
   *<outFilePrefix_ann1>*, *<outFilePrefix_ann2>* the names to be prefixed to all the output files generated for annotaion1 and annotation2 respectively.  
  
    See help() for more details, optional parameters, their defaults and output files generated.

**3. Integration and p-value calculation :** Examine enrichment based on the value assigned to the 
calibration parameter alpha (indicates the relative weight of the SNP annotation1). Call *integrationPvalCalc()* using the tansformed zscores obtained in the previous step. This function returns 
the p-value calculated for the alpha. 

   Enrichment can also be examined by varying the relative importance of the two annotation datasets, 
i.e. for a list of alpha values, by calling *integrationPvalCalc()* multiple times e.g with lapply()

   #### Example
   *pvalue <- integrationPvalCalc(<rdaAnn1>, <rdaAnn2>, <outFile_prefix>, <alpha>)*  
   where  
   *<rdaAnn1>* and *<rdaAnn2>* are the transformed z-score rda file names for annotation1 and annotation2 obtained from the Randomimzation and Transformation step.
   *<outFile_prefix>* the name to be prefixed to all the output files generated.  
   *<alpha>* is the relative weight to be assigned to the annotation1 can have values ranging from 0 to 1.default is 0.5  
   
   See help() for details on further optional parameters, their defaults and output files generated.
