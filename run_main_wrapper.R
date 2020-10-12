source("Rscripts/coloc_project_functions.R")

main(eqtl="data/eQTL_CAU_Height.txt.gz",gwas="data/GWAS_CAU_Height.txt.gz",
     mode="bse",
     eqtlBetaCol=5,eqtlSeCol=6,
     gwasBetaCol=2,gwasSeCol=3)