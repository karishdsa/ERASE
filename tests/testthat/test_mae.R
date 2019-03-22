context("ERASE-MAE")

test_that("ERASE-MAE works", {
  load("testDataASE.rda")
  load("testDataGWAS.rda")
  load("testDataAnnSpi.rda")

  #Intersection
  asegwas <- getIntersection(dfAseData, dfGwasData)
  expect_true(is.data.frame(asegwas))
  asegwas$neglogpval <- -log10(asegwas$p)


  aseann <- getIntersection(dfAseData,
                            dfAnnData,
                            snpAnnType="SPI",
                            snpAnnCol2Sel = "abs_dpsi_max_tissue")
  expect_true(is.data.frame(aseann))


  ase_gwas_ann <- getIntersectionMae(asegwas, aseann[ , c("cmp.col", "abs_dpsi_max_tissue")])
  expect_true(is.data.frame(ase_gwas_ann))

  #Randomization & transformation
  asegwas_zscoreFile <- randomization(df_sigASE_SNPann = ase_gwas_ann[ ase_gwas_ann$min.FDR < 0.05, ],
                    df_nonASE_SNPann = ase_gwas_ann[ ase_gwas_ann$min.FDR >= 0.05, ],
                    colname_rankSNPann = "neglogpval",
                    colname_chk4distr = "avgReads",
                    outFilePrefix= "MAE_AseGwas",
                    nIterations = 5 ,
                    binwidth = 20,
                    mode="MAE",
                    seedValue=12345 )

  aseann_zscoreFile <- randomization(df_sigASE_SNPann = ase_gwas_ann[ ase_gwas_ann$min.FDR < 0.05, ],
                                      df_nonASE_SNPann = ase_gwas_ann[ ase_gwas_ann$min.FDR >= 0.05, ],
                                      colname_rankSNPann = "abs_dpsi_max_tissue",
                                      colname_chk4distr = "avgReads",
                                      outFilePrefix= "MAE_AseAnn",
                                      nIterations = 5 ,
                                      binwidth = 20,
                                      mode="MAE",
                                      seedValue=67890 )

    #Integration and p-value calculation
  pval <- integrationPvalCalc(rdaSnpAnn1 = asegwas_zscoreFile,
                                rdaSnpAnn2 = aseann_zscoreFile,
                                outFilePrefix = "Inte_AseGwasAnn",
                                alphaVal = 0.5,
                                seedValue = 1234)

})
