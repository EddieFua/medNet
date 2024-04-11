rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Second data/clean data.RData')
setwd('/Users/fuyinghao/Documents/Mediator analysis/save data')
library(medNet)
name = c('Trans fats','Calories from fat','Calories from sweets, desserts')
sel_exposure = c('FF_DT_TRFAT','FF_PCTFAT','FF_PCTSWEET')
sel_exposure = toupper(sel_exposure)
idx = match(sel_exposure, colnames(exposure))
name = c('Trans fats','Calories from fat','Calories from sweets, desserts')
edgelist = data.frame(NA,NA,NA)
colnames(edgelist) = c("from" ,"to" ,  "type")
for (j in idx){
  load(paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/',colnames(mediator)[j],'_original_res.RData'))
  final_net = original_res$final_net
  final_net[,1] = c(rep(name[which(idx==j)],length(original_res$final_net[original_res$final_net[,1]=='<NA>',1])),final_net[!is.na(original_res$final_net[,1]),1])
  edgelist = rbind(edgelist,final_net)
}
edgelist = edgelist[complete.cases(edgelist),]
write.table(edgelist,file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/all_metabolites.txt',row.names = F,sep = '/',quote = F)

i = 3
edgelist = data.frame(NA,NA,NA)
colnames(edgelist) = c("from" ,"to" ,  "type")
for (j in idx){
  load( paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_',i,'_',j,'.RData'))
  final_net = disturbed_res$final_net
  final_net[,1] = c(rep(name[which(idx==j)],length(disturbed_res$final_net[disturbed_res$final_net[,1]=='<NA>',1])),final_net[!is.na(disturbed_res$final_net[,1]),1])
  edgelist = rbind(edgelist,final_net)
}
edgelist = edgelist[complete.cases(edgelist),]
write.table(edgelist,file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_all_metabolites.txt',row.names = F,sep = '/',quote = F)



pathway = data.frame(NA,NA,NA,NA,NA,NA)
colnames(pathway) = c("pathway",        "num_in_selected", "num_in_all",      "p_value",         "compound id",'mediator'    )
for (j in idx){
  load(paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/',colnames(mediator)[j],'_original_res.RData'))
  pathway = rbind(pathway,cbind(original_pathway,mediator = rep(name[which(idx==j)],nrow(original_pathway))))
}
pathway = pathway[complete.cases(pathway),]
write.csv(pathway,'/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/orginal_pathway.csv')


pathway_disturbed = data.frame(NA,NA,NA,NA,NA,NA)
colnames(pathway_disturbed) = c("pathway",        "num_in_selected", "num_in_all",      "p_value",         "compound id",'mediator'    )
for (j in idx){
  load( paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_',3,'_',j,'.RData'))
  pathway_disturbed = rbind(pathway_disturbed,cbind(disturbed_pathway,mediator = rep(name[which(idx==j)],nrow(disturbed_pathway))))
}
pathway_disturbed = pathway_disturbed[complete.cases(pathway_disturbed),]
write.csv(pathway_disturbed,'/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/disturbed_pathway.csv')








