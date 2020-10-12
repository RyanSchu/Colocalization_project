source("Rscripts/coloc_project_functions.R")

main(eqtl="data/eQTL_CAU_Height.txt.gz",gwas="data/GWAS_CAU_Height.txt.gz",
     mode="bse",
     eqtlGeneCol=1,eqtlSNPCol=2,eqtlBetaCol=5,eqtlSeCol=6,eqtlSampleSize=578,
     gwasSNPCol=1,gwasBetaCol=2,gwasSeCol=3,gwasSampleSize=49796,gwasMAFCol=4,
     outFile="test_CAU_out.txt")