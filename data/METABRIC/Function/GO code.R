rm(list = ls())
library(medNetmy)
load('/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/cleandata.RData')
load('/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/res.RData')
library("org.Hs.eg.db")
library(GOstats)
library(Category)
library(qvalue)
res = res_p
res$final_net[,1] = ifelse(is.na(res$final_net[,1]=='<NA>'),1,res$final_net[,1])
mediator = c(res$final_net[,1],res$final_net[,2])[-(1:106)]
info <- AnnotationDbi::select(org.Hs.eg.db, keys=mediator,
                              columns=c("ENSEMBL"), 
                              keytype="SYMBOL")
ens <- info$ENSEMBL
str_name = strsplit(colnames(dat$mediator),split = '|',fixed = T)
for (i in 1:ncol(dat$mediator)) {
  colnames(dat$mediator)[i] = str_name[[i]][1]
}
names = AnnotationDbi::select(org.Hs.eg.db, keys=colnames(dat$mediator),
                              columns=c("ENSEMBL"), 
                              keytype="SYMBOL")
df<-data.frame(name = names$ENSEMBL)

analyze<- function(ens,df){
  sel.entrez<-  unique(na.omit(ens)) # selected gene names (15 character)
  all.entrez<-  unique(na.omit(df$name)) # all gene names in the data matrix under study
  
  # if there is more than one matched, then all are taken to the next step
  
  sel.entrez<-unlist(mget(substr(sel.entrez,1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))
  sel.entrez<-unique(sel.entrez[!is.na(sel.entrez)])
  all.entrez<-unlist(mget(substr(all.entrez,1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))
  all.entrez<-unique(all.entrez[!is.na(all.entrez)])
  
  # this is hypergeometric testing, by GOstats
  
  params <- new("GOHyperGParams", geneIds=sel.entrez[!is.na(sel.entrez)], universeGeneIds=all.entrez[!is.na(all.entrez)], ontology="BP", pvalueCutoff=0.05,conditional=T, testDirection="over", annotation="org.Hs.eg.db")
  
  over = hyperGTest(params)
  ov<-summary(over)
  ov<-ov[ov[,6]<=500 & ov[,6]>=12,]
  
  for(i in 2:4) ov[,i]<-signif(ov[,i], 3)
  
  # this is the fill in which genes caused the pathway to be significant           
  ov<-cbind(ov, idx = rep('',nrow(ov)), name = rep('', nrow(ov)))
  ov[,8]<-as.vector(ov[,8])
  ov[,9]<-as.vector(ov[,9])
  
  for(i in 1:nrow(ov))
  {
    this.genes<-mget(ov[i,1],org.Hs.egGO2ALLEGS)[[1]]
    this.genes<-unique(this.genes[this.genes %in% sel.entrez])
    that<-this.genes[1]
    if(length(this.genes)>1)
    {
      for(j in 2:length(this.genes)) that<-paste(that, ", ", this.genes[j], sep="")
    }
    ov[i,8]<-that
    
    this.genes<-unlist(mget(this.genes, org.Hs.egSYMBOL))
    that<-this.genes[1]
    if(length(this.genes)>1)
    {
      for(j in 2:length(this.genes)) that<-paste(that, ", ", this.genes[j], sep="")
    }
    ov[i,9]<-that
    
  }
  ov <- cbind(ov, p_adjusted = rep(0,nrow(ov)))
  ov$p_adjusted[which(is.na(ov$Pvalue)==FALSE)] <- p.adjust(ov$Pvalue,"BH")
  return(ov)
  
}
res <- analyze(ens,df)
ov = res
save(ov, file = '/Users/fuyinghao/Documents/Mediator analysis/metabric res/pathway.RData')
write.csv(ov,file = '/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/GO term_p.csv')


