rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/clean data.RData')
setwd('/Users/fuyinghao/Documents/Mediator analysis/save data')
library(medNet)
name = c('Trans fats','Calories from fat','Calories from sweets, desserts')
sel_exposure = c('FF_DT_TRFAT','FF_PCTFAT','FF_PCTSWEET')
sel_exposure = toupper(sel_exposure)
idx = match(sel_exposure, colnames(exposure))
enrichment_analysis = function(mediator, res){
  meta_name_sub = unique(colnames(mediator))
  meta_name_100 = unique(c(res$final_net[,2],na.omit(res$final_net[,1])))
  # meta_name_100 = unique(c(res_p$final_net[,2],res_p$final_net[c(31:35,63),1]))
  library(metapone)
  data(pa)
  pathway_name_100 <- c()
  
  for (i in 1:length(meta_name_100)) {
    p <- pa[which(pa$KEGG.ID==meta_name_100[i]),2]
    pathway_name_100 <- c(pathway_name_100, p)
  }
  
  # pathway related to the top 100 compound
  pathway_name_100 = unique(pathway_name_100)
  
  df = as.data.frame(matrix(nrow=0,ncol=5))
  for (n in pathway_name_100) {
    c = pa[which(pa$pathway.name==n),4]
    in_100 = sum(c %in% meta_name_100)
    in_sub = sum(c %in% meta_name_sub)
    in_100_name = paste(c[c %in% meta_name_100],collapse=",")
    in_sub_name = paste(c[c %in% meta_name_sub],collapse=",")
    # p_value = 1.0-phyper(in_100, length(imp), length(meta_name_sub)-length(imp), in_sub)
    p_value = 1.0-phyper(in_100-1, in_sub, length(meta_name_sub)-in_sub, length(meta_name_100))
    df <- rbind(df, c(n, in_100_name, in_sub_name, p_value, paste(c[c %in% meta_name_100],collapse=",")))
  }
  
  names(df) <- c("pathway",  "num_in_selected", "num_in_all", "p_value", "compound id")
  sum(df$p_value<0.01)
  df = df[df$p_value<0.01,]
  df = df[order(df$p_value),]
  return(df)
}
idx1_matrix = matrix(NA,ncol = 4,nrow = 10)
idx2_matrix = matrix(NA,ncol = 4,nrow = 10)
idx3_matrix = matrix(NA,ncol = 4,nrow = 10)
idx1_matrix_meta = matrix(NA,ncol = 4,nrow = 10)
idx2_matrix_meta = matrix(NA,ncol = 4,nrow = 10)
idx3_matrix_meta = matrix(NA,ncol = 4,nrow = 10)
for (m in 1:10){
for(i in idx){
  original_res = findGlobalNet(response = response,mediator = as.matrix(mediator),exposure = as.matrix(exposure[,i]), mediator.network = mediator_network,H=10,T = 0)
  original_pathway = enrichment_analysis(mediator,original_res)
  save(original_pathway,original_res,file = paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/',colnames(mediator)[i],'_original_res.RData'))
}

disturbed_mediator_network = mediator_network
sample_idx = c()
num_samples = floor(0.1 * nrow(mediator_network))
new_idx = sample(seq_len(nrow(mediator_network)), size = num_samples)
sample_idx = c(sample_idx,new_idx)
for (i in seq_along(new_idx)) {
  disturb_edge = disturbed_mediator_network[new_idx[i], ]
  available_edges = which((!(disturbed_mediator_network[, 2] %in% disturb_edge))&(!(disturbed_mediator_network[, 1] %in% disturb_edge)))
  available_edges = available_edges[which(!available_edges%in%sample_idx)]
  exchange_idx = sample(available_edges, 1)
  sample_idx = c(sample_idx,exchange_idx)
  temp = disturbed_mediator_network[exchange_idx, 2]
  disturbed_mediator_network[exchange_idx, 2] = disturbed_mediator_network[new_idx[i], 2]
  disturbed_mediator_network[new_idx[i], 2] = temp
}
k = 1
for(j in idx){
  disturbed_res = findGlobalNet(response = response,mediator = as.matrix(mediator),exposure = as.matrix(exposure[,j]), mediator.network = disturbed_mediator_network,H=10,T = 0)
  disturbed_pathway = enrichment_analysis(mediator,disturbed_res)
  save(disturbed_res,disturbed_pathway,file = paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_',k,'_',j,'.RData'))
}

for(k in 2:4) {
  id = (1:nrow(disturbed_mediator_network))[!(1:nrow(disturbed_mediator_network))  %in% sample_idx]
  new_idx = sample(id, size = num_samples)
  sample_idx = c(sample_idx, new_idx)
  for (i in seq_along(new_idx)) {
    disturb_edge = disturbed_mediator_network[new_idx[i], ]
    available_edges = which((!(disturbed_mediator_network[, 2] %in% disturb_edge))&(!(disturbed_mediator_network[, 1] %in% disturb_edge)))
    available_edges = available_edges[which(!available_edges%in%sample_idx)]
    exchange_idx = sample(available_edges, 1)
    sample_idx = c(sample_idx,exchange_idx)
    temp = disturbed_mediator_network[exchange_idx, 2]
    disturbed_mediator_network[exchange_idx, 2] = disturbed_mediator_network[new_idx[i], 2]
    disturbed_mediator_network[new_idx[i], 2] = temp
  }
  for (j in idx) {
    disturbed_res = findGlobalNet(
      response = response,
      mediator = as.matrix(mediator),
      exposure = as.matrix(exposure[, j]),
      mediator.network = disturbed_mediator_network,
      H = 10,
      T = 0
    )
    disturbed_pathway = enrichment_analysis(mediator, disturbed_res)
    save(
      disturbed_res,
      disturbed_pathway,
      file = paste0(
        '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_',
        k,
        '_',
        j,
        '.RData'
      )
    )
  }
}
# save(new_idx_matrix,file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_idx.RData')
intersect_pathway_num = matrix(NA,nrow = length(idx),ncol = 4)
intersect_meta_num = matrix(NA,nrow = length(idx),ncol = 4)
for(i in 1:4){
  for (j in idx){
    load( paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_',i,'_',j,'.RData'))
    load(paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/',colnames(mediator)[j],'_original_res.RData'))
    intersect_pathway_num[which(idx==j),i] = length(intersect(disturbed_pathway$pathway,original_pathway$pathway))
    intersect_meta_num[which(idx==j),i] = length(intersect(unique(disturbed_res$final_net[,2]), unique(original_res$final_net[,2])))
  }}
idx1_matrix[m,] = intersect_pathway_num[1,]
idx2_matrix[m,] = intersect_pathway_num[2,]
idx3_matrix[m,] = intersect_pathway_num[3,]
idx1_matrix_meta[m,] = intersect_meta_num[1,]
idx2_matrix_meta[m,] = intersect_meta_num[2,]
idx3_matrix_meta[m,] =intersect_meta_num[3,]
}

colMeans(idx1_matrix,na.rm = T)
colMeans(idx2_matrix,na.rm = T)
colMeans(idx3_matrix,na.rm = T)
colMeans(idx1_matrix_meta,na.rm = T)
colMeans(idx2_matrix_meta,na.rm = T)
colMeans(idx3_matrix_meta,na.rm = T)
save(idx1_matrix,idx2_matrix,idx3_matrix,idx1_matrix_meta,idx2_matrix_meta,idx3_matrix_meta,file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/ten_times.RData')



load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/ten_times.RData')
