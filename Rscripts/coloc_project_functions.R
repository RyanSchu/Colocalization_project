suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(dplyr))

"%&%" = function(a,b) paste(a,b,sep="")

check.fread<-function(file,header=T,stringsAsFactors=F){ 
  #Author Ryan Schubert
  #Meant to be run in unix environments with the zcat program installed
  #Not for the R environ, but useful for the client
  if (grepl("\\.gz$",file)){
    fread("zcat " %&% file,header=header,stringsAsFactors = stringsAsFactors)
  } else {
    fread(file,header=header,stringsAsFactors = stringsAsFactors)
  }
}

get_gene_list<-function(eqtl_df,column=1){ 
  #Author Ryan Schubert
  #takes in the data frame of eqtl summary stats and 
  #extracts the column with the gene names in them, 
  #returns as a vector
  eqtl_df[,column] %>% unlist %>% unname %>% unique()
}

get_gene_eqtl_effects<-function(gene_id,eqtl_df,betaCol=5,seCol=6,geneCol=1){
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the effect sizes and variances for each snp
  #If the gene is not found in the eqtl df then return an empty list
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id,]
  if (nrow(eqtl_df) == 0) { return(list())}
  beta<-eqtl_df[,betaCol] %>% unlist %>% unname
  var<-eqtl_df[,seCol]^2 %>% unlist %>% unname
  return(list(beta=beta,var=var))
}

get_gene_eqtl_pvalue<-function(gene_id,eqtl_df,pvalCol=4,geneCol=1){ 
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the pvalue of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id,]
  if (nrow(eqtl_df) == 0) { return(list())} 
  p<-eqtl_df[,pvalCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gene_eqtl_maf<-function(gene_id,eqtl_df,mafCol=3,geneCol=1){ 
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the maf of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the eqtl data, though it is common practice
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id,]
  if (nrow(eqtl_df) == 0) { return(list())} 
  p<-eqtl_df[,mafCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gwas_effects<-function(gwas_df,betaCol=5,seCol=6,geneCol=1){
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the effect sizes and variances for each snp
  #If the gene is not found in the gwas df then return an empty list
  beta<-gwas_df[,betaCol] %>% unlist %>% unname
  var<-gwas_df[,seCol]^2 %>% unlist %>% unname
  return(list(beta=beta,var=var))
}

get_gwas_pvalue<-function(gwas_df,pvalCol=4,geneCol=1){ 
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the pvalue of each snp
  #If the gene is not found in the gwas df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  p<-gwas_df[,pvalCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gwas_maf<-function(gwas_df,mafCol=3,geneCol=1){ 
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the maf of each snp
  #If the gene is not found in the gwas df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the gwas data, though it is common practice
  p<-gwas_df[,mafCol] %>% unlist %>% unname
  return(list(p=p))
}

main<-function(eqtl,gwas,mode="bse",
               eqtlGeneCol=NULL,eqtlSNPCol=NULL,eqtlMAFCol=NULL,eqtlPvalCol=NULL,eqtlBetaCol=NULL,eqtlSeCol=NULL,
               gwasSNPCol=NULL,gwasMAFCol=NULL,gwasPvalCol=NULL,gwasBetaCol=NULL,gwasSeCol=NULL){
  if( mode  == "bse") {
    cat("Running coloc using the beta/standard error method\n")
    if(is.null(eqtlBetaCol) && is.null(eqtlSeCol)) {print("eQTL Beta/Se column not set. Please set these and try again. Exiting"); return() }
    if(is.null(gwasBetaCol) && is.null(gwasSeCol)) {print("GWAS Beta/Se column not set. Please set these and try again. Exiting"); return() }
  } else if ( mode == "p" ){
    cat("Running coloc using the pvalue method\n")
    if(is.null(eqtlPvalCol)) {print("eQTL Pval column not set. Please set this and try again. Exiting"); return() }
    if(is.null(gwasPvalCol)) {print("GWAS Pval column not set. Please set this and try again. Exiting"); return() }
  } else {
    cat("Mode not recognized. Please choose one of [ bse, p ]. Exiting\n"); return()
  }
  
  #uncomment the follwing lines when running in a Unix/linux environ
  # cat("Reading in data\n")
  # eqtldf<-check.fread(eqtl,header=T)
  # gwasdf<-check.fread(gwas,header=T)
  
  #comment out the following lines when running in Unix/Linux environ
  cat("Reading in data\n")
  eqtldf<-read.table(eqtl,header=T)
  gwasdf<-read.table(gwas,header=T)
  
  cat("Getting gene list\n")
  gene_list<-get_gene_list(eqtldf)
  ngenes<-length(gene_list)
  if ( ngenes == 0 ) {print("0 genes found. Exiting.\n"); return()}
  cat(ngenes,"genes found\n")
  
  if ( mode == "bse") {
    QTL_Betas<-get_gene_eqtl_effects
    
  }
  
}

