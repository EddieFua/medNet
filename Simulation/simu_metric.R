rm(list = ls())
setwd('/Users/fuyinghao/Documents/Mediator analysis/stimu res')
source('/Users/fuyinghao/Documents/Mediator analysis/StimulateData.R')
library(dplyr)
library(medNetmy)
combn = c(200,400,600)
result_matrix_2 = matrix(NA,ncol = 6,nrow = length(combn))
for (k in 1:length(combn)) {
  e_TP = c()
  e_FP = c()
  m_TP = c()
  m_FP = c()
  edge_TP = c()
  edge_FP = c()
  for (j in 1:10){
    n_exposure = 5
    n_mediator = 900
    true_exposure = 1
    data = StimulateData_update(
      n_exposure = n_exposure,
      n_mediator = n_mediator,
      N = combn[k],
      ba_m = 3,
      ba_power = 0.5,
      multi_rate = 0.5,
      poisson_rate = 1,
      max_mediator_link = 10,
      true_exposure = true_exposure,
      num_link = 10)
    res = medNetmy::findGlobalNet(data$response,data$exposure,data$mediator,data$mediator_network,T = 0.6,H = 14,method = 'LR')
    
    ##exposure 
    res_exp = unique(res$final_net[res$final_net[,3] == 'mediation_link',1])
    res_exp = substr(res_exp,2,nchar(res_exp))
    # TP_exp = sum(res_exp%in%unique(data$link$mediation_link[,1]))/(sum(res_exp%in%unique(data$link$mediation_link[,1]))+sum(!unique(data$link$mediation_link[,1])%in%res_exp))
    TP_exp = sum(res_exp%in%unique(data$link$mediation_link[,1]))/true_exposure
    # FP_exp = sum(!res_exp%in%unique(data$link$mediation_link[,1]))/sum(!(1:n_exposure)%in%unique(data$link$mediation_link[,1]))
    FP_exp = sum(!res_exp%in%unique(data$link$mediation_link[,1]))/(ncol(data$exposure)-true_exposure)
    
    ##mediator TPR
    ##include those mediators in the functional link
    res_med = unique(c(res$final_net[,2],res$final_net[res$final_net$type=='functional_link',1],res$final_net[res$final_net$type=='functional_link',2]))
    # sum(!res$final_net[res$final_net$type=='functional_link',1]%in%res$final_net[,2])
    ##exclude those mediators in the functional link
    res_med = substr(res_med,2,nchar(res_med))
    true_med = unique(which(data$exposure_mediator==1,arr.ind = T)[,2])
    TP_mediator = sum(res_med%in%true_med)/length(true_med)
    FP_mediator = sum(!res_med%in%true_med)/sum(!(1:n_mediator)%in%true_med)

    
    ###edge
    true_edge = which(data$exposure_mediator==1,arr.ind = T)
    # true_edge = data$link$mediation_link
    med_link = res$final_net[res$final_net$type=='mediation_link',1:2]
    med_link = t(apply(med_link, 1,function(i)as.numeric(substr(i,2,nchar(i)))))
    mm_name = substr(res$final_net[res$final_net$type=='mediation_link',2],2,nchar(res$final_net[res$final_net$type=='mediation_link',2]))
    m_name_1 = substr(res$final_net[res$final_net$type=='functional_link',2],2,nchar(res$final_net[res$final_net$type=='functional_link',2]))
    m_name_2 = substr(res$final_net[res$final_net$type=='functional_link',1],2,nchar(res$final_net[res$final_net$type=='functional_link',1]))
    for (i in 1:nrow(med_link)) {
      id_1 = substr(res$final_net[res$final_net$type=='functional_link',1],2, nchar(res$final_net[res$final_net$type=='functional_link',1]))%in%med_link[i,2]
      id_2 = substr(res$final_net[res$final_net$type=='functional_link',2],2, nchar(res$final_net[res$final_net$type=='functional_link',2]))%in%med_link[i,2]
      if(sum(id_1)>0|sum(id_2)>0){
        med_link = rbind(med_link,matrix(c(rep(med_link[i,1],sum(id_1)),m_name_1[id_1]),ncol = 2))
        med_link = rbind(med_link,matrix(c(rep(med_link[i,1],sum(id_2)),m_name_2[id_2]),ncol = 2))
      }
    }
    med_link = unique(apply(med_link, 2, as.numeric))
    TP_edge = sum(duplicated(rbind(med_link,true_edge)))/nrow(true_edge)
    FP_edge = (nrow(med_link)-sum(duplicated(rbind(med_link,true_edge))))/(n_exposure*n_mediator-nrow(true_edge))
    e_TP = c(e_TP,TP_exp)
    e_FP = c(e_FP,FP_exp)
    m_TP = c(m_TP,TP_mediator)
    m_FP = c(m_FP,FP_mediator)
    edge_TP = c(edge_TP,TP_edge)
    edge_FP = c(edge_FP,FP_edge)
  }
  result_matrix_2[k,] = c(mean(e_TP),mean(e_FP),mean(m_TP),mean(m_FP),mean(edge_TP),mean(edge_FP))
}

save(result_matrix_2,file= '/Users/fuyinghao/Documents/Mediator analysis/stimu res/result_matrix.RData')








