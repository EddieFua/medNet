rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/clean data.RData')
setwd('/Users/fuyinghao/Documents/Mediator analysis/save data')
library(medNet)
disturb_network <- function(network, propor) {
  num_samples <- floor(propor * nrow(network))
  sample_idx <- sample(seq_len(nrow(network)), size = num_samples)
  
  for (i in seq_along(sample_idx)) {
    # Identify the edge to be disturbed
    disturb_edge <- network[sample_idx[i], ]
    
    # Use vectorized comparison to avoid apply()
    available_edges <- which((!(network[-sample_idx, 2] %in% disturb_edge))&(!(network[-sample_idx, 1] %in% disturb_edge)))
    
    exchange_idx <- sample(available_edges, 1)
    
    # Swap the edges
    temp <- network[exchange_idx, 2]
    network[exchange_idx, 2] <- network[sample_idx[i], 2]
    network[sample_idx[i], 2] <- temp
  }
  
  return(network)
}
name = c('Trans fats','Calories from fat','Calories from sweets, desserts')
sel_exposure = c('FF_DT_TRFAT','FF_PCTFAT','FF_PCTSWEET')
# sel_exposure = c('FF_DT_POTA','FF_DT_CALC','FF_DT_malt','FF_free_choline','FF_phosphocholine','FF_DT_lact')
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

for(i in idx){
  original_res = findGlobalNet(response = response,mediator = as.matrix(mediator),exposure = as.matrix(exposure[,i]), mediator.network = mediator_network,H=10,T = 0)
  original_pathway = enrichment_analysis(mediator,original_res)
  save(original_pathway,original_res,file = paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/',colnames(mediator)[i],'_original_res.RData'))
}


disturbed_mediator_network = mediator_network
sample_idx = c()
num_samples = floor(0.05 * nrow(mediator_network))
new_idx = sample(seq_len(nrow(mediator_network)), size = num_samples)
sample_idx = c(sample_idx,new_idx)
new_idx_matrix = matrix(ncol = num_samples,nrow = 12)
new_idx_matrix[1,] = new_idx
for (i in seq_along(new_idx)) {
  disturb_edge = disturbed_mediator_network[new_idx[i], ]
  available_edges = which((!(disturbed_mediator_network[-sample_idx, 2] %in% disturb_edge))&(!(disturbed_mediator_network[-sample_idx, 1] %in% disturb_edge)))
  exchange_idx = sample(available_edges, 1)
  new_idx_matrix[2,i] = exchange_idx
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

for(k in 2:6) {
  id = (1:nrow(mediator_network))[!(1:nrow(mediator_network))  %in% sample_idx]
  new_idx = sample(id, size = num_samples)
  new_idx_matrix[2*k-1,] = new_idx 
  sample_idx = c(sample_idx, new_idx)
  for (i in seq_along(new_idx)) {
    disturb_edge = disturbed_mediator_network[new_idx[i],]
    available_edges = which((!(
      disturbed_mediator_network[-new_idx, 2] %in% disturb_edge
    )) & (!(
      disturbed_mediator_network[-new_idx, 1] %in% disturb_edge
    )))
    exchange_idx = sample(available_edges, 1)
    new_idx_matrix[2*k,i] = exchange_idx
    sample_idx = c(sample_idx, exchange_idx)
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
save(new_idx_matrix,file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_idx.RData')
intersect_pathway_num = matrix(NA,nrow = length(idx),ncol = 6)
intersect_meta_num = matrix(NA,nrow = length(idx),ncol = 6)
for(i in 1:6){
for (j in idx){
  load( paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_',i,'_',j,'.RData'))
  load(paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/',colnames(mediator)[j],'_original_res.RData'))
  intersect_pathway_num[which(idx==j),i] = length(intersect(disturbed_pathway$pathway,original_pathway$pathway))
  intersect_meta_num[which(idx==j),i] = length(intersect(unique(disturbed_res$final_net[,2]), unique(original_res$final_net[,2])))
}}

disturbed_rate = seq(0.05,0.3,0.05)*2

data = data.frame(disturbed_rate, intersect_meta_num = intersect_meta_num[1,])

# metarate the line plot for metas
meta_plot_1 <- ggplot(data, aes(x = disturbed_rate, y = intersect_meta_num)) +
  geom_line(color = 'blue') + # Line color
  geom_point(color = 'red') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection meta Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Metabolites")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


# Assuming 'intersect_pathway_num' is another vector corresponding to 'disturbed_rate'
# Create a data frame for pathway data
pathway_data <- data.frame(disturbed_rate, intersect_pathway_num = intersect_pathway_num[1,])

# metarate the line plot for pathways
pathway_plot_1 <- ggplot(pathway_data, aes(x = disturbed_rate, y = intersect_pathway_num)) +
  geom_line(color = 'green') + # Line color
  geom_point(color = 'orange') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection Pathway Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Pathways")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


data = data.frame(disturbed_rate, intersect_meta_num = intersect_meta_num[2,])

# metarate the line plot for metas
meta_plot_2 <- ggplot(data, aes(x = disturbed_rate, y = intersect_meta_num)) +
  geom_line(color = 'blue') + # Line color
  geom_point(color = 'red') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection meta Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Metabolites")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


# Assuming 'intersect_pathway_num' is another vector corresponding to 'disturbed_rate'
# Create a data frame for pathway data
pathway_data <- data.frame(disturbed_rate, intersect_pathway_num = intersect_pathway_num[2,])

# metarate the line plot for pathways
pathway_plot_2 <- ggplot(pathway_data, aes(x = disturbed_rate, y = intersect_pathway_num)) +
  geom_line(color = 'green') + # Line color
  geom_point(color = 'orange') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection Pathway Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Pathways")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


data = data.frame(disturbed_rate, intersect_meta_num = intersect_meta_num[3,])

# metarate the line plot for metas
meta_plot_3 <- ggplot(data, aes(x = disturbed_rate, y = intersect_meta_num)) +
  geom_line(color = 'blue') + # Line color
  geom_point(color = 'red') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection meta Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Metabolites")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


# Assuming 'intersect_pathway_num' is another vector corresponding to 'disturbed_rate'
# Create a data frame for pathway data
pathway_data <- data.frame(disturbed_rate, intersect_pathway_num = intersect_pathway_num[3,])

# metarate the line plot for pathways
pathway_plot_3 <- ggplot(pathway_data, aes(x = disturbed_rate, y = intersect_pathway_num)) +
  geom_line(color = 'green') + # Line color
  geom_point(color = 'orange') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection Pathway Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Pathways")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )

cowplot::plot_grid(pathway_plot_1,meta_plot_1,pathway_plot_2,meta_plot_2,pathway_plot_3,meta_plot_3,ncol = 2,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/disturbed_network_meta.png', width = 15, height = 18)

