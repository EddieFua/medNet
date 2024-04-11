rm(list = ls())
setwd('/Users/fuyinghao/Documents/Mediator analysis/save data')
library(medNet)
load('/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/cleandata.RData')
intersect_num = matrix(NA,ncol = 11,nrow = 10)
intersect_gene_num = matrix(NA,ncol = 11,nrow = 10)
ov_less_num = matrix(NA,ncol = 11,nrow = 10)
ov_num = matrix(NA,ncol = 11,nrow = 10)

for (z in seq(400, 550, 50)) {
  intersect_tmp = c()
  ov_less_tmp = c()
  ov_tmp =  c()
  intersect_gene_tmp = c()
  for (j in 1:10){
  idx_0 = sample(which(dat$response == 0), z)
  idx_1 = sample(which(dat$response == 1), z)
  new_dat = list()
  new_dat$response = as.matrix(c(dat$response[idx_0], dat$response[idx_1]), ncol = 1)
  new_dat$exposure = as.matrix(c(dat$exposure[idx_0], dat$exposure[idx_1]), ncol = 1)
  new_dat$mediator = rbind(dat$mediator[idx_0, ], dat$mediator[idx_1, ])
  new_dat$mediator_network = dat$mediator_network
  res = medNet::findGlobalNet(
    response = new_dat$response,
    mediator = new_dat$mediator,
    exposure = new_dat$exposure,
    mediator.network = new_dat$mediator_network,
    H = 14,
    T = 0.5
  )
  
  library("org.Hs.eg.db")
  library(GOstats)
  library(Category)
  library(qvalue)
  res$final_net[, 1] = ifelse(is.na(res$final_net[, 1] == '<NA>'), 1, res$final_net[, 1])
  mediator = unique(res$final_net[, 2])
  info <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = mediator,
    columns = c("ENSEMBL"),
    keytype = "SYMBOL"
  )
  ens <- info$ENSEMBL
  str_name = strsplit(colnames(dat$mediator), split = '|', fixed = T)
  for (i in 1:ncol(dat$mediator)) {
    colnames(dat$mediator)[i] = str_name[[i]][1]
  }
  names = AnnotationDbi::select(
    org.Hs.eg.db,
    keys = colnames(dat$mediator),
    columns = c("ENSEMBL"),
    keytype = "SYMBOL"
  )
  df <- data.frame(name = names$ENSEMBL)
  
  analyze <- function(ens, df) {
    sel.entrez <-
      unique(na.omit(ens)) # selected gene names (15 character)
    all.entrez <-
      unique(na.omit(df$name)) # all gene names in the data matrix under study
    
    # if there is more than one matched, then all are taken to the next step
    
    sel.entrez <-
      unlist(mget(
        substr(sel.entrez, 1, 15),
        org.Hs.egENSEMBL2EG,
        ifnotfound = NA
      ))
    sel.entrez <- unique(sel.entrez[!is.na(sel.entrez)])
    all.entrez <-
      unlist(mget(
        substr(all.entrez, 1, 15),
        org.Hs.egENSEMBL2EG,
        ifnotfound = NA
      ))
    all.entrez <- unique(all.entrez[!is.na(all.entrez)])
    
    # this is hypergeometric testing, by GOstats
    
    params <-
      new(
        "GOHyperGParams",
        geneIds = sel.entrez[!is.na(sel.entrez)],
        universeGeneIds = all.entrez[!is.na(all.entrez)],
        ontology = "BP",
        pvalueCutoff = 0.05,
        conditional = T,
        testDirection = "over",
        annotation = "org.Hs.eg.db"
      )
    
    over = hyperGTest(params)
    ov <- summary(over)
    ov <- ov[ov[, 6] <= 500 & ov[, 6] >= 12, ]
    
    for (i in 2:4)
      ov[, i] <- signif(ov[, i], 3)
    
    # this is the fill in which genes caused the pathway to be significant
    ov <- cbind(ov, idx = rep('', nrow(ov)), name = rep('', nrow(ov)))
    ov[, 8] <- as.vector(ov[, 8])
    ov[, 9] <- as.vector(ov[, 9])
    
    for (i in 1:nrow(ov))
    {
      this.genes <- mget(ov[i, 1], org.Hs.egGO2ALLEGS)[[1]]
      this.genes <- unique(this.genes[this.genes %in% sel.entrez])
      that <- this.genes[1]
      if (length(this.genes) > 1)
      {
        for (j in 2:length(this.genes))
          that <- paste(that, ", ", this.genes[j], sep = "")
      }
      ov[i, 8] <- that
      
      this.genes <- unlist(mget(this.genes, org.Hs.egSYMBOL))
      that <- this.genes[1]
      if (length(this.genes) > 1)
      {
        for (j in 2:length(this.genes))
          that <- paste(that, ", ", this.genes[j], sep = "")
      }
      ov[i, 9] <- that
      
    }
    ov <- cbind(ov, p_adjusted = rep(0, nrow(ov)))
    ov$p_adjusted[which(is.na(ov$Pvalue) == FALSE)] <-
      p.adjust(ov$Pvalue, "BH")
    return(ov)
    
  }
  res_pathway <- analyze(ens, df)
  res_pathway = res_pathway[res_pathway$Pvalue<0.03,]
  ov_less = res_pathway
  res_less = res
  load(
    '/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/pathway_update.RData'
  )
  load(
    '/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/res_update.RData'
  )
  intersect_tmp = c(intersect_tmp,length(intersect(ov_less$Term, ov$Term)))
  ov_less_tmp = c(ov_less_tmp,sum(!ov_less$Term%in%intersect(ov_less$Term, ov$Term)))
  ov_tmp =  c(ov_tmp,sum(!ov$Term%in%intersect(ov_less$Term, ov$Term)))
  intersect_gene_tmp = c(intersect_gene_tmp,length(intersect(unique(res_less$final_net[,2]), unique(res$final_net[,2]))))
  }
  intersect_num[,which(z==seq(50, 550, 50))] = intersect_tmp
  ov_less_num[,which(z==seq(50, 550, 50))] = ov_less_tmp
  ov_num[,which(z==seq(50, 550, 50))] =  ov_tmp
  intersect_gene_num[,which(z==seq(50, 550, 50))] = intersect_gene_tmp
}


save(intersect_num,ov_less_num,ov_num,intersect_gene_num,file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/less_sample.RData')
library(ggplot2)


intersect_long <- reshape2::melt(intersect_num)
mapping <- setNames(as.character(seq(100, 1100, 100)), as.character(1:11))
intersect_long$Var2 <- factor(mapping[intersect_long$Var2], levels = as.character(seq(100, 1100, 100)))
p1 = ggplot(intersect_long, aes(x = as.factor(Var2), y = value)) +
  geom_boxplot(aes(fill = as.factor(Var2)), show.legend = FALSE) +
  stat_summary(fun.y = mean, geom = "point", aes(group = 1),
               shape = 18, size = 3, color = "red") +
  stat_summary(fun.y = mean, geom = "line", aes(group = 1),
               color = "blue", size = 1) +
  labs(title = "Boxplot Analysis with Connected Means:",
       subtitle = 'Intersected Pathway Terms Across Varying Sample Sizes',
       x = "Sample Size",
       y = "Number of Intersect Terms") +
  theme_minimal() +
  theme_bw()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))+
  scale_x_discrete(breaks = seq(100, 1100, 100), labels = seq(100, 1100, 100))

intersect_gene_long <- reshape2::melt(intersect_gene_num)
mapping <- setNames(as.character(seq(100, 1100, 100)), as.character(1:11))
intersect_gene_long$Var2 <- factor(mapping[intersect_gene_long$Var2], levels = as.character(seq(100, 1100, 100)))
p2 = ggplot(intersect_gene_long, aes(x = as.factor(Var2), y = value)) +
  geom_boxplot(aes(fill = as.factor(Var2)), show.legend = FALSE) +
  stat_summary(fun.y = mean, geom = "point", aes(group = 1),
               shape = 18, size = 3, color = "red") +
  stat_summary(fun.y = mean, geom = "line", aes(group = 1),
               color = "blue", size = 1) +
  labs(title = "Boxplot Analysis with Connected Means:",
       subtitle = 'Intersected Gene Counts Across Varying Sample Sizes',
       x = "Sample Size",
       y = "Number of Intersect Genes") +
  theme_minimal() +
  theme_bw()+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))+
  scale_x_discrete(breaks = seq(100, 1100, 100), labels = seq(100, 1100, 100))
cowplot::plot_grid(p1,p2,ncol = 2,nrow = 1)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/less_sample.png', width = 10, height = 5)






library(VennDiagram)
fill_colors <- c("#56B4E9", "#E69F00")
border_colors <- c("black", "black")
cat_col <- c("darkblue", "darkorange")



# 生成Venn图
venn.plot <- venn.diagram(
  x = list(Pathway1 = unique(res_less$final_net[,2]), Pathway2 = unique(res$final_net[,2])),
  category.names = c("Generated With Less Sample", "Generated with All Sample"),
  filename = "/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/less sample_gene.png",
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
  main.pos = c(0.5, 0.9),  # 调整标题位置
  lwd = 2,
  # height = 5,
  # width = 10,
  resolution = 300
)

# 生成Venn图
venn.plot <- venn.diagram(
  x = list(Pathway1 = ov_less$Term, Pathway2 = ov$Term),
  category.names = c("Generated With Less Sample", "Generated with All sample"),
  filename = "/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/less sample network.png",
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








