rm(list = ls())
years = 5
dat = as.data.frame(read.csv('/Users/fuyinghao/Documents/Msc T1/ML for biomedical research/Assignment/Assignment1/data.csv'))
dat[dat=='']=NA
###delete samples who died of other causes
dat = dat[dat$death_from_cancer!="Died of Other Causes",]
dat = dat[-which(is.na(dat$death_from_cancer)),]
###censor
dat$overall_survival_months = dat$overall_survival_months/12
dat = dat[-which(dat[which(dat$death_from_cancer=="Living"),]$overall_survival_months<years),]
write.csv(dat,file = '/Users/fuyinghao/Documents/Mediator analysis/metabric res/METRABIC_clean.csv')
clinical_data = dat[,1:31]
geneexp_data = dat[,32:520]
###get mediator network
response = clinical_data$overall_survival
exposure = as.matrix(clinical_data$nottingham_prognostic_index,ncol = 1)
mediator = geneexp_data
gene_link = read.table('/Users/fuyinghao/Documents/Mediator analysis/LUAD data/HomoSapiens_binary_hq.txt',na.strings = "NA",sep = '\t',header = T)
Rownames_med = sapply(1:length(colnames(mediator)),function(x) stringr::str_to_upper(colnames(mediator)[x]),simplify = 'vector')
colnames(mediator) = Rownames_med
mediator_network_1 = match(gene_link$Gene_A,Rownames_med)
mediator_network_2 = match(gene_link$Gene_B,Rownames_med)
mediator_network = as.data.frame(cbind(mediator_network_1,mediator_network_2))
mediator_network = na.omit(mediator_network)
dat = list(response = response, exposure= exposure, mediator = mediator,mediator_network = mediator_network)
save(dat,file = '/Users/fuyinghao/Documents/Mediator analysis/metabric res/cleandata.RData')

rm(list = ls())
library(medNetmy)
load('/Users/fuyinghao/Documents/Mediator analysis/metabric res/cleandata.RData')
setwd('/Users/fuyinghao/Documents/Mediator analysis/metabric res')
s_p = Sys.time()
res_p = medNetmy::findGlobalNet(response = dat$response,mediator = dat$mediator,exposure = dat$exposure, mediator.network = dat$mediator_network,
                    H=18,T = 0.6)
e_p = Sys.time()
e_p-s_p
res = res_p
save(res,file = '/Users/fuyinghao/Documents/Mediator analysis/metabric res/res.RData')
s = Sys.time()
mediator_network = dat$mediator_network
res = medNet::findGlobalNet(response = dat$response,mediator = dat$mediator,exposure = dat$exposure, mediator.network = dat$mediator_network,
                              H=17,T = 0.6)
e = Sys.time()
res$AUC







