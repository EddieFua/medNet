rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/clean data.RData')
setwd('/Users/fuyinghao/Documents/Mediator analysis/save data')
library(medNetmy)
# exposure = exposure[,na.omit(match(unique(res_p$final_net[,1]),colnames(exposure)))[-1]]
exposure = exposure[,-8]
res = list()
for (i in 1:ncol(exposure)) {
  e = as.matrix(exposure[,i])
  colnames(e) = colnames(exposure)[i]
  err = try(medNetmy::findGlobalNet(response = as.integer(response),mediator = mediator,exposure = e, mediator.network = mediator_network,
                          H=14 ,T = 0.59,method = 'LR'),TRUE)
  if('try-error'%in% class(err))
  {
    err <- NA
  }
  res[[i]] <- err
}

id = which(sapply(1:length(res),function(i)is.list(res[[i]])))
res_p = list()
for (i in 1:length(id)) {
  res_p[[i]] = res[[id[i]]]
}
res = res_p
# res = sapply(1:ncol(exposure), function(i)medNetmy::findGlobalNet(response = as.integer(response),mediator = mediator,exposure = as.matrix(exposure[,i]), mediator.network = mediator_network,
# H=12 ,T = 0.6))

# load(file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/res.RData')
pathway_frame = matrix(NA,ncol = 4)[-1,]
for (j in 1:length(id)) {
  meta_name_sub = unique(colnames(mediator))
  meta_name_100 = unique(c(res[[j]]$final_net[,2],na.omit(res[[j]]$final_net[,1])))
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
  pathway_frame = rbind(pathway_frame,cbind(as.matrix(rep(colnames(exposure)[id[j]],nrow(df))),df[,c(1,2,4)]))
}

# save(pathway_frame,file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/df.RData')
# write.csv(pathway_frame,file = '/Users/fuyinghao/Documents/Mediator analysis/Second data/pathway_frame.csv')



load('/Users/fuyinghao/Documents/Mediator analysis/Second data/df.RData')
colnames(pathway_frame) = c('exposure',colnames(pathway_frame)[-1])
edge_list = cbind(pathway_frame$exposure,pathway_frame$pathway)
tmp = graph_from_edgelist(edge_list)
tmp = as.undirected(tmp)
type = c()
for (i in V(tmp)$name) {
  if(i %in%colnames(exposure)){
    type = c(type,1)
  }else{
    type = c(type,2)
  }
}
V(tmp)$type = type
V(tmp)$color <- c("Orange", "gray70")[V(tmp)$type]
V(tmp)$name[which(V(tmp)$type == 1)] = c('Calcium','Potassium','Lactose','Maltose','Free choline','Phosphocholine')
V(tmp)$shape <- c("circle","vrectangle")[V(tmp)$type]
V(tmp)$size <- ifelse(V(tmp)$type == 1, 26, 22)
V(tmp)$label.cex <- ifelse(igraph::V(tmp)$type == 1,0.85,0.8) 
layout = layout.fruchterman.reingold(tmp)
# V(tmp)$name = gsub(' ','\n',V(tmp)$name)
e = get.edgelist(tmp,names=FALSE)
l = qgraph::qgraph.layout.fruchtermanreingold(e,vcount=vcount(tmp),
                                              area=20*(vcount(tmp)^2),repulse.rad=(vcount(tmp)^3.1),niter = 2000)
pdf('/Users/fuyinghao/Documents/Mediator analysis/Second data/All_pathway_glm.pdf',height = 10,width = 10)
plot(tmp,layout =  l,width = 1,edge.color = 'black',vertex.color = adjustcolor("SkyBlue2", alpha.f = .5),vertex.label.color= "black")
dev.off()

library(ggraph)
library(tidygraph)
pdf('/Users/fuyinghao/Documents/Mediator analysis/Second data/pathway_ggraph.pdf',width = 20,height = 10)
tmp%>%
  ggraph(layout = 'linear', circular = TRUE) +
  # ggraph('lgl')+
  geom_edge_arc(colour= "gray50",
                lineend = "round",
                strength = .1) +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 point.padding = unit(0.2, "lines"),size = 5,
                 colour="gray10") +
  geom_node_point(color = V(tmp)$color, size = 5) +
  scale_edge_width(range = c(0, 2.5)) +
  scale_edge_alpha(range = c(0, .2)) +
  theme_graph(background = "white") +
  theme(legend.position = "top") +
  guides(edge_width = FALSE,
         edge_alpha = FALSE)
dev.off()





