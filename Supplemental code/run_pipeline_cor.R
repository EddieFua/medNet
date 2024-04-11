run_pipeline_cor <- function(H,path,num_link,method) {
  source('/Users/fuyinghao/Documents/Mediator analysis/Revision/code/simu_cor.R')
  library(dplyr)
  library(medNet)
  combn = c(200, 300, 400)
  result_matrix_2 = matrix(NA, ncol = 5, nrow = 54)
  T_c = H
  H_c = c(14, 14, 14)
  e_TP_matrix = matrix(nrow = 54, ncol = 10)
  e_FP_matrix = matrix(nrow = 54, ncol = 10)
  m_TP_matrix = matrix(nrow = 54, ncol = 10)
  m_FP_matrix = matrix(nrow = 54, ncol = 10)
  m_FDR_matrix = matrix(nrow = 54, ncol = 10)
  edge_TP_matrix = matrix(nrow = 54, ncol = 10)
  edge_FP_matrix = matrix(nrow = 54, ncol = 10)
  edge_FDR_matrix = matrix(nrow = 54, ncol = 10)
  num_mediator_matrix = matrix(nrow = 54, ncol = 10)
  for (j in 1:10) {
    n_exposure = 3
    n_mediator = 500
    true_exposure = 1
    data = StimulateData_update(
      n_exposure = n_exposure,
      n_mediator = n_mediator,
      N = combn,
      ba_m = 3,
      ba_power = 0.5,
      multi_rate = 0.8,
      poisson_rate = 1,
      max_mediator_link = 10,
      true_exposure = true_exposure,
      num_link = num_link
    )
    save(data,file = paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/simu_data/',j,'.RData'))
    e_TP = c()
    e_FP = c()
    m_TP = c()
    m_FP = c()
    m_FDR = c()
    edge_TP = c()
    edge_FP = c()
    edge_FDR = c()
    num_mediator = c()
    for (m in method){
      for (q in 1:length(T_c)){
        for (k in 1:6) {
          res = medNet::findGlobalNet(
            data[[k]]$response,
            data[[k]]$exposure,
            data[[k]]$mediator,
            data[[k]]$mediator_network,
            T = T_c[q],
            H = 14,
            method = m,
            fold.number = 10,
            thread = 4
          )
          
          ##exposure
          res_exp = unique(res$final_net[res$final_net[, 3] == 'mediation_link', 1])
          res_exp = substr(res_exp, 2, nchar(res_exp))
          TP_exp = sum(res_exp %in% unique(data[[k]]$link$mediation_link[, 1])) / true_exposure
          FP_exp = sum(!res_exp %in% unique(data[[k]]$link$mediation_link[, 1])) /
            (ncol(data[[k]]$exposure) - true_exposure)
          
          ##mediator TPR
          ##include those mediators in the functional link
          res_med = unique(c(res$final_net[, 2], res$final_net[res$final_net$type ==
                                                                 'functional_link', 1], res$final_net[res$final_net$type == 'functional_link', 2]))
          ##exclude those mediators in the functional link
          res_med = substr(res_med, 2, nchar(res_med))
          true_med = unique(which(data[[k]]$exposure_mediator == 1, arr.ind = T)[, 2])
          TP_mediator = sum(res_med %in% true_med) / length(true_med)
          FP_mediator = sum(!res_med %in% true_med) / sum(!(1:n_mediator) %in% true_med)
          
          
          ###edge
          true_edge = which(data[[k]]$exposure_mediator == 1, arr.ind = T)
          # true_edge = data$link$mediation_link
          med_link = res$final_net[res$final_net$type == 'mediation_link', 1:2]
          med_link = t(apply(med_link, 1, function(i)
            as.numeric(substr(i, 2, nchar(
              i
            )))))
          mm_name = substr(res$final_net[res$final_net$type == 'mediation_link', 2], 2, nchar(res$final_net[res$final_net$type ==
                                                                                                              'mediation_link', 2]))
          m_name_1 = substr(res$final_net[res$final_net$type == 'functional_link', 2], 2, nchar(res$final_net[res$final_net$type ==
                                                                                                                'functional_link', 2]))
          m_name_2 = substr(res$final_net[res$final_net$type == 'functional_link', 1], 2, nchar(res$final_net[res$final_net$type ==
                                                                                                                'functional_link', 1]))
          for (i in 1:nrow(med_link)) {
            id_1 = substr(res$final_net[res$final_net$type == 'functional_link', 1], 2, nchar(res$final_net[res$final_net$type ==
                                                                                                              'functional_link', 1])) %in% med_link[i, 2]
            id_2 = substr(res$final_net[res$final_net$type == 'functional_link', 2], 2, nchar(res$final_net[res$final_net$type ==
                                                                                                              'functional_link', 2])) %in% med_link[i, 2]
            if (sum(id_1) > 0 | sum(id_2) > 0) {
              med_link = rbind(med_link, matrix(c(
                rep(med_link[i, 1], sum(id_1)), m_name_1[id_1]
              ), ncol = 2))
              med_link = rbind(med_link, matrix(c(
                rep(med_link[i, 1], sum(id_2)), m_name_2[id_2]
              ), ncol = 2))
            }
          }
          med_link = unique(apply(med_link, 2, as.numeric))
          TP_edge = sum(duplicated(rbind(med_link, true_edge))) / nrow(true_edge)
          FP_edge = (nrow(med_link) - sum(duplicated(rbind(
            med_link, true_edge
          )))) / (n_exposure * n_mediator - nrow(true_edge))
          
          
          ###mediator FDR
          FDR = sum(!res_med %in% true_med) / length(res_med)
          e_FDR = (nrow(med_link) - sum(duplicated(rbind(
            med_link, true_edge
          )))) / nrow(med_link)
          num_mediator = c(num_mediator, sum(data[[k]]$exposure_mediator == 1))
          
          e_TP = c(e_TP, TP_exp)
          e_FP = c(e_FP, FP_exp)
          m_TP = c(m_TP, TP_mediator)
          m_FP = c(m_FP, FP_mediator)
          m_FDR = c(m_FDR, FDR)
          edge_TP = c(edge_TP, TP_edge)
          edge_FP = c(edge_FP, FP_edge)
          edge_FDR = c(edge_FDR, e_FDR)
          
        }}}
    e_TP_matrix[, j] = e_TP
    e_FP_matrix[, j] = e_FP
    m_TP_matrix[, j] = m_TP
    m_FP_matrix[, j] = m_FP
    m_FDR_matrix[, j] = m_FDR
    edge_TP_matrix[, j] = edge_TP
    edge_FP_matrix[, j] = edge_FP
    edge_FDR_matrix[,j ] = edge_FDR
    num_mediator_matrix[, j] = num_mediator
  }
  result_matrix_2[,1 ] = rowMeans(e_TP_matrix)
  result_matrix_2[,2 ] = rowMeans(e_FP_matrix)
  result_matrix_2[,3 ] = rowMeans(m_TP_matrix)
  result_matrix_2[,4 ] = rowMeans(m_FP_matrix)
  result_matrix_2[,5 ] = rowMeans(m_FDR_matrix)
  save(
    num_mediator_matrix,
    e_TP_matrix,
    e_FP_matrix,
    m_TP_matrix,
    m_FP_matrix,
    m_FDR_matrix,
    edge_TP_matrix,
    edge_FP_matrix,
    edge_FDR_matrix,
    file = paste0(
      path,
      'scenario2.RData'
    )
  )
  return(result_matrix_2)
}
