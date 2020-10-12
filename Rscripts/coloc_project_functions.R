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


get_gene_list<-function(eqtl_df,geneCol=1){ 
  #Author Ryan Schubert
  #takes in the data frame of eqtl summary stats and 
  #extracts the column with the gene names in them, 
  #returns as a vector
  eqtl_df[,geneCol] %>% unlist %>% unname %>% unique()
}


get_eqtl_snps<-function(gene_id,eqtl_df,snpCol=2,geneCol=1){ 
  #Author Ryan Schubert
  #take in the eqtl data frame
  #returns the snps in the data frame
  #If the gene is not found in the eqtl df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the eqtl data, though it is common practice
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id,]
  if (nrow(eqtl_df) == 0) { return(list())} 
  snp<-eqtl_df[,snpCol] %>% unlist %>% unname
  return(list(snp=snp))
}

get_gene_eqtl_effects<-function(gene_id,eqtl_df,betaCol=5,seCol=6,geneCol=1,snpCol,snpList){
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the effect sizes and variances for each snp
  #If the gene is not found in the eqtl df then return an empty list
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id & eqtl_df[,snpCol] %in% snpList,]
  if (nrow(eqtl_df) == 0) { return(list())}
  beta<-eqtl_df[,betaCol] %>% unlist %>% unname
  var<-eqtl_df[,seCol]^2 %>% unlist %>% unname
  return(list(beta=beta,var=var))
}

get_gene_eqtl_pvalue<-function(gene_id,eqtl_df,pvalCol=4,geneCol=1,snpCol,snpList){ 
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the pvalue of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id & eqtl_df[,snpCol] %in% snpList,]
  if (nrow(eqtl_df) == 0) { return(list())} 
  p<-eqtl_df[,pvalCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gene_eqtl_maf<-function(gene_id,eqtl_df,mafCol=3,geneCol=1,snpCol,snpList){ 
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the maf of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the eqtl data, though it is common practice
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id & eqtl_df[,snpCol] %in% snpList,]
  if (nrow(eqtl_df) == 0) { return(list())} 
  maf<-eqtl_df[,mafCol] %>% unlist %>% unname
  return(list(maf=maf))
}


get_gwas_snps<-function(gwas_df,snpCol=1){ 
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the snps in the data frame
  #If the gene is not found in the gwas df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the gwas data, though it is common practice
  snp<-gwas_df[,snpCol] %>% unlist %>% unname
  return(list(snp=snp))
}

get_gwas_effects<-function(gwas_df,betaCol,seCol,snpCol,snpList){
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the effect sizes and variances for each snp
  #If the gene is not found in the gwas df then return an empty list
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if (nrow(gwas_df) == 0) { return(list())} 
  beta<-gwas_df[,betaCol] %>% unlist %>% unname
  var<-gwas_df[,seCol]^2 %>% unlist %>% unname
  return(list(beta=beta,var=var))
}

get_gwas_pvalue<-function(gwas_df,pvalCol,snpCol,snpList){ 
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the pvalue of each snp
  #If the gene is not found in the gwas df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if (nrow(gwas_df) == 0) { return(list())} 
  p<-gwas_df[,pvalCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gwas_maf<-function(gwas_df,mafCol,snpCol,snpList){ 
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the maf of each snp
  #If the gene is not found in the gwas df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the gwas data, though it is common practice
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if (nrow(gwas_df) == 0) { return(list())} 
  maf<-gwas_df[,mafCol] %>% unlist %>% unname
  return(list(maf=maf))
}



main<-function(eqtl,gwas,mode="bse",
               eqtlGeneCol=NULL,eqtlSNPCol=NULL,eqtlMAFCol=NULL,eqtlPvalCol=NULL,eqtlBetaCol=NULL,eqtlSeCol=NULL, eqtlSampleSize=NULL,
               gwasSNPCol=NULL,gwasMAFCol=NULL,gwasPvalCol=NULL,gwasBetaCol=NULL,gwasSeCol=NULL,gwasSampleSize=NULL,
               outFile=NULL){
  if(is.null(eqtlSampleSize) || is.null(gwasSampleSize)) {print("Sample sizes not set. Please make sure both GWAS and QTL sample sizes are set. Exiting."); return()}
  if( mode  == "bse") {
    cat("Running coloc using the beta/standard error method\n")
    if(is.null(eqtlBetaCol) || is.null(eqtlSeCol)) {print("eQTL Beta/Se column not set. Please set these and try again. Exiting"); return() }
    if(is.null(gwasBetaCol) || is.null(gwasSeCol)) {print("GWAS Beta/Se column not set. Please set these and try again. Exiting"); return() }
  } else if ( mode == "p" ){
    cat("Running coloc using the pvalue method\n")
    if(is.null(eqtlPvalCol)) {print("eQTL Pval column not set. Please set this and try again. Exiting"); return() }
    if(is.null(gwasPvalCol)) {print("GWAS Pval column not set. Please set this and try again. Exiting"); return() }
  } else {
    cat("Mode not recognized. Please choose one of [ bse, p ]. Exiting\n"); return()
  }
  
  #uncomment the follwing lines when running in a Unix/linux environ
  # cat("Reading in data\n")
  # eqtldf<-check.fread(eqtl,header=T,stringsAsFactors = F)
  # gwasdf<-check.fread(gwas,header=T,stringsAsFactors = F)
  
  #comment out the following lines when running in Unix/Linux environ
  cat("Reading in data\n")
  eqtldf<-read.table(eqtl,header=T,stringsAsFactors = F)
  gwasdf<-read.table(gwas,header=T,stringsAsFactors = F)
  
  cat("Getting gene list\n")
  gene_list<-get_gene_list(eqtldf)
  ngenes<-length(gene_list)
  if ( ngenes == 0 ) {print("0 genes found. Exiting.\n"); return()}
  cat(ngenes,"genes found\n")
  if ( !is.null(outFile) ){
    names<-list(nsnps=numeric(),PP.H0.abf=numeric(),PP.H1.abf=numeric(),PP.H2.abf=numeric(),PP.H3.abf=numeric(),PP.H4.abf=numeric(),gene=character())
    write.table(names,outFile,sep='\t',quote=F)
  }
  
  if ( mode == "bse") {
    cat("Processing GWAS\n")
     
    GWAS_snps<-get_gwas_snps(gwas_df = gwasdf,snpCol=gwasSNPCol)
  
    for (i in 1:ngenes){
      gene<-gene_list[i]
      cat("Processing gene",gene," ",i,"/",ngenes,"\n")
      QTL_snps<-get_eqtl_snps(gene_id=gene,eqtl_df=eqtldf,snpCol=eqtlSNPCol,geneCol=eqtlGeneCol)
      intersection<-base::intersect(GWAS_snps$snp,QTL_snps$snp)
      nsnps<-length(intersection)
      if (nsnps == 0) { print("0 snps present in intersection for this gene. skipping"); next}
      cat(nsnps, "snps found in GWAS and eQTL intersection for this gene\n")
      GWAS_effects<-get_gwas_effects(gwas_df = gwasdf,betaCol = gwasBetaCol,seCol=gwasSeCol,snpCol=gwasSNPCol,snpList=intersection)
      QTL_effects<-get_gene_eqtl_effects(gene_id=gene,eqtl_df = eqtldf,betaCol=eqtlBetaCol,seCol=eqtlSeCol,snpCol=eqtlSNPCol,snpList=intersection)
      if (length(GWAS_effects$beta) != length(QTL_effects$beta)) { print("List of effect sizes of differing lengths. Duplicates SNPs may be present. Please resolve and rerun. Exiting."); return()}
      
      if (!is.null(gwasMAFCol)){ 
        print("using maf from gwas data set")
        maf<-get_gwas_maf(gwas_df=gwasdf,mafCol=gwasMAFCol,snpCol = gwasSNPCol,snpList=intersection)
      } else if(!is.null(eqtlMAFCol)){
        print("using maf from QTL data set")
        maf<-get_gene_eqtl_maf(gene_id=gene,eqtl_df = eqtldf,mafcol=eqtlMAFCol,geneCol=eqtlGeneCol,snpCol=eqtlSNPCol,snpList=intersection)
      }
      coloc_result <- coloc.abf(dataset1=list(beta=GWAS_effects$beta, varbeta=GWAS_effects$var, N=gwasSampleSize,type="quant"),
                          dataset2=list(beta=QTL_effects$beta, varbeta=QTL_effects$var, N=eqtlSampleSize,type="quant"),
                          MAF=maf$maf)
      summary<-coloc_result$summary
      summary$gene<-gene
      summary<-bind_cols(summary)
      str(summary)
      if ( !is.null(outFile) ){
        write.table(summary,outFile,append=T,sep='\t',col.names=F,quote=F,row.names=F)
      }
#      print(summary)
    }
  } else if ( mode == "p"){
    print("Not done yet")
    return()
  }
  
}
?write.table
