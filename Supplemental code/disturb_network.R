rm(list = ls())
library(medNet)
load('/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/cleandata.RData')
disturb_edge = sample(1:nrow(dat$mediator_network),floor(nrow(dat$mediator_network)*0.1))
node = unique(c(dat$mediator_network[,1],dat$mediator_network[,2]))
for (i in disturb_edge){
  sample_node = sample(node,1)
  if (sample_node != dat$mediator_network[i,1] ){
    dat$mediator_network[i,2] = sample_node
  }
}
dat$mediator_network = dat$mediator_network[-disturb_edge,]

setwd('/Users/fuyinghao/Documents/Mediator analysis/save data')
res = findGlobalNet(response = dat$response,mediator = dat$mediator,exposure = dat$exposure, mediator.network = dat$mediator_network,
                                H=18,T = 0.6)

library("org.Hs.eg.db")
library(GOstats)
library(Category)
library(qvalue)
res$final_net[,1] = ifelse(is.na(res$final_net[,1]=='<NA>'),1,res$final_net[,1])
mediator = unique(res$final_net[,2])
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
res_pathway <- analyze(ens,df)
ov = res_pathway
save(ov,res, file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/pathway_normal.RData')

ov_disturbed = ov
# load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/pathway_delete.RData')


rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/pathway_deleted.RData')
ov_deleted = ov
res_deleted = res
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/pathway_disturbed.RData')
ov_disturbed = ov
res_disturbed = res
# load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/pathway_normal.RData')
load('/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/pathway_update.RData')

library(VennDiagram)
fill_colors <- c("#56B4E9", "#E69F00")
border_colors <- c("black", "black")
cat_col <- c("darkblue", "darkorange")

# 生成Venn图
venn.plot <- venn.diagram(
  x = list(Pathway1 = unique(res_disturbed$final_net[,2]), Pathway2 = unique(res$final_net[,2])),
  category.names = c("Generated With Disturbed Network", "Generated with Normal Network"),
  filename = "/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/disturbed_gene.png",
  output = TRUE,
  imagetype = "png",
  fill = fill_colors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "serif",
  cat.cex = 1.5,
  cat.col = cat_col,
  cat.fontfamily = "serif",
  cat.pos = c(0, 180),  # 调整标签位置
  cat.dist = 0.02,      # 调整标签距离圆心的距离
  margin = 0.05,  # 减少边距
  main = "Venn Diagram of Genes",
  main.cex = 1.5,  # 减小标题大小
  main.fontfamily = "serif",
  main.col = "black",
  main.just = "center",  # 标题居中对齐
  main.pos = c(0.5, 1),  # 调整标题位置
  lwd = 2,
  # height = 5,
  # width = 10,
  resolution = 300
)

venn.plot <- venn.diagram(
  x = list(Pathway1 = unique(res_deleted$final_net[,2]), Pathway2 = unique(res$final_net[,2])),
  category.names = c("Generated With Deleted Network", "Generated with Normal Network"),
  filename = "/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/deleted_gene.png",
  output = TRUE,
  imagetype = "png",
  fill = fill_colors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "serif",
  cat.cex = 1.5,
  cat.col = cat_col,
  cat.fontfamily = "serif",
  cat.pos = c(0, 180),  # 调整标签位置
  cat.dist = 0.02,      # 调整标签距离圆心的距离
  margin = 0.05,  # 减少边距
  main = "Venn Diagram of Genes",
  main.cex = 1.5,  # 减小标题大小
  main.fontfamily = "serif",
  main.col = "black",
  main.just = "center",  # 标题居中对齐
  main.pos = c(0.5, 1),  # 调整标题位置
  lwd = 2,
  # height = 5,
  # width = 10,
  resolution = 300
)

# 生成Venn图
venn.plot <- venn.diagram(
  x = list(Pathway1 = ov_disturbed$Term, Pathway2 = ov$Term),
  category.names = c("Generated by Deleted Network", "Generated by Original Network"),
  filename = "/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/Deleted network.png",
  output = TRUE,
  imagetype = "png",
  fill = fill_colors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "serif",
  cat.cex = 1.5,
  cat.col = cat_col,
  cat.fontfamily = "serif",
  cat.pos = c(0, 180),  # 调整标签位置
  cat.dist = 0.02,      # 调整标签距离圆心的距离
  margin = 0.05,  # 减少边距
  main = "Venn Diagram of Pathway Terms",
  main.cex = 1.5,  # 减小标题大小
  main.fontfamily = "serif",
  main.col = "black",
  main.just = "center",  # 标题居中对齐
  main.pos = c(0.5, 0.9),  # 调整标题位置
  lwd = 2,
  # height = 5,
  # width = 10,
  resolution = 300
)


# 生成Venn图
venn.plot <- venn.diagram(
  x = list(Pathway1 = ov_disturbed$Term, Pathway2 = ov$Term),
  category.names = c("Generated by Disturbed Network", "Generated by Original Network"),
  filename = "/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/Disturbed network.png",
  output = TRUE,
  imagetype = "png",
  fill = fill_colors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "serif",
  cat.cex = 1.5,
  cat.col = cat_col,
  cat.fontfamily = "serif",
  cat.pos = c(0, 180),  # 调整标签位置
  cat.dist = 0.02,      # 调整标签距离圆心的距离
  margin = 0.05,  # 减少边距
  main = "Venn Diagram of Pathway Terms",
  main.cex = 1.5,  # 减小标题大小
  main.fontfamily = "serif",
  main.col = "black",
  main.just = "center",  # 标题居中对齐
  main.pos = c(0.5, 0.9),  # 调整标题位置
  lwd = 2,
  # height = 5,
  # width = 10,
  resolution = 300
)




