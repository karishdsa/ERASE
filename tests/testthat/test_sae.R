context("ERASE-SAE")

test_that("ERASE-SAE works", {
  load("testDataASE.rda")
  load("testDataGWAS.rda")
  asegwas <- getIntersection(dfAseData, dfGwasData)
  expect_true(is.data.frame(asegwas))

  asegwas$neglogpval <- -log10(asegwas$p)
  pval <- randomization(df_sigASE_SNPann = asegwas[ asegwas$min.FDR < 0.05, ],
                df_nonASE_SNPann = asegwas[ asegwas$min.FDR >= 0.05, ],
                colname_rankSNPann = "neglogpval",
                colname_chk4distr = "avgReads",
                outFilePrefix= "SAE",
                nIterations = 5 ,
                binwidth = 10,
                mode="SAE",
                seedValue=12345 )

  expect_equal(pval, 0.02299062)
})
