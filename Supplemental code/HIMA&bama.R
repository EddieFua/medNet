rm(list = ls())
library(hdmed)
library(HIMA)
library(bama)
library(medNet)
setwd('/Users/fuyinghao/Documents/Mediator analysis/stimu res')
HIMA_TPR_matrix = matrix(nrow = 20, ncol = 9)
HIMA_FPR_matrix = matrix(nrow = 20, ncol = 9)
HIMA_FDR_matrix = matrix(nrow = 20, ncol = 9)
medfix_TPR_matrix = matrix(nrow = 20, ncol = 9)
medfix_FPR_matrix = matrix(nrow = 20, ncol = 9)
medfix_FDR_matrix = matrix(nrow = 20, ncol = 9)
bama_TPR_matrix = matrix(nrow = 20, ncol = 9)
bama_FPR_matrix = matrix(nrow = 20, ncol = 9)
bama_FDR_matrix = matrix(nrow = 20, ncol = 9)
mednet_TPR_matrix = matrix(nrow = 20, ncol = 9)
mednet_FPR_matrix = matrix(nrow = 20, ncol = 9)
mednet_FDR_matrix = matrix(nrow = 20, ncol = 9)
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/benchmark.RData')
H_v = c(16,16,18,16,16,18,16,16,18)
for (j in 1:20) {
  print(j)
  n_exposure = 3
  n_mediator = 500
  true_exposure = 1
  load(paste0('/Users/fuyinghao/Documents/Mediator analysis/Revision/simu_data/nonlinear_',j,'.RData'))
  HIMA_TPR = c()
  HIMA_FPR = c()
  HIMA_FDR = c()
  medfix_TPR = c()
  medfix_FPR = c()
  medfix_FDR = c()
  bama_TPR = c()
  bama_FPR = c()
  bama_FDR = c()
  mednet_TPR = c()
  mednet_FPR = c()
  mednet_FDR = c()
    for (k in 1:6) {
        
      
        p <- ceiling(nrow(data[[k]]$mediator) / log(nrow(data[[k]]$mediator)))
        true_med = unique(which(data[[k]]$exposure_mediator == 1, arr.ind = T)[, 2])
        # #HIMA
        # out_hima <-
        #   hima(X = as.vector(data[[k]]$exposure[,1]),
        #     Y = as.vector(data[[k]]$response), as.matrix(data[[k]]$mediator),as.matrix(data[[k]]$exposure[,2:3]),
        #     topN = p,
        #     Y.family = "gaussian",
        #     M.family = "gaussian",
        #     penalty = "MCP",
        #     max.iter = 10^5,
        #     parallel = T
        #   )
        # #HIMA metric
        # hima_mediator = rownames(out_hima)
        # hima_mediator = as.numeric(gsub("`", "", hima_mediator))
        # # hima_mediator = as.numeric(gsub("\\D", "", hima_mediator))
        # HIMA_TPR = c(HIMA_TPR,sum(hima_mediator %in% true_med) / length(true_med))
        # HIMA_FPR = c(HIMA_FPR,sum(!hima_mediator %in% true_med) / sum(!(1:n_mediator) %in% true_med))
        # HIMA_FDR = c(HIMA_FDR,sum(!hima_mediator %in% true_med) / length(hima_mediator))
        
        
      #medfix
        out_medfix = mediate_medfix(
          as.vector(data[[k]]$exposure), as.matrix(data[[k]]$mediator),
          as.vector(data[[k]]$response),
          C1 = as.matrix(data[[k]]$exposure[,1]), C2 = as.matrix(data[[k]]$mediator),
          nlambda = 100,
          nlambda2 = 50,
          nfolds = 10,
          seed = 1
        )
        medfix_mediator = (1:500)[out_medfix$contributions$beta_pv<=0.05]
        medfix_TPR = c(medfix_TPR,sum(medfix_mediator %in% true_med) / length(true_med))
        medfix_FPR = c(medfix_FPR,sum(!medfix_mediator %in% true_med) / sum(!(1:n_mediator) %in% true_med))
        medfix_FDR = c(medfix_FDR,sum(!medfix_mediator %in% true_med) / length(medfix_mediator))
        
      #bama
      #   out_bama = bama(as.vector(data[[k]]$response),as.vector(data[[k]]$exposure[,1]), as.matrix(data[[k]]$mediator),C1 = as.matrix(data[[k]]$exposure[,1]), C2 = as.matrix(data[[k]]$mediator),method = "BSLMM",burnin = 160, ndraws = 175,
      #     control = list(k = 2, lm0 = 1e-04, lm1 = 1, lma1 = 1, l = 1))
      # out_bama = summary(out_bama)
      # bama_mediator = (1:500)[out_bama[,4]!=0]
      # bama_TPR = c(bama_TPR,sum(bama_mediator %in% true_med) / length(true_med))
      # bama_FPR = c(bama_FPR,sum(!bama_mediator %in% true_med) / sum(!(1:n_mediator) %in% true_med))
      # bama_FDR = c(bama_FDR,sum(!bama_mediator %in% true_med) / length(bama_mediator))

      # #mednet
      # res = medNet::findGlobalNet(
      #   data[[k]]$response,
      #   as.matrix(data[[k]]$exposure[,1]),
      #   data[[k]]$mediator,
      #   data[[k]]$mediator_network,
      #   T = 0.6,
      #   H = H_v[k],
      #   method = 'LR',
      #   fold.number = 10,
      #   thread = 2
      # )
      # res_med = unique(c(res$final_net[, 2], res$final_net[res$final_net$type ==
      #                                                        'functional_link', 1], res$final_net[res$final_net$type == 'functional_link', 2]))
      # res_med = substr(res_med, 2, nchar(res_med))
      # mednet_TPR = c(mednet_TPR,sum(res_med %in% true_med) / length(true_med))
      # mednet_FPR = c(mednet_FPR,sum(!res_med %in% true_med) / sum(!(1:n_mediator) %in% true_med))
      # mednet_FDR = c(mednet_FDR,sum(!res_med %in% true_med) / length(res_med))
    }
  # HIMA_TPR_matrix[j,] =HIMA_TPR
  # HIMA_FPR_matrix[j,] =HIMA_FPR
  # HIMA_FDR_matrix[j,] =HIMA_FDR
  medfix_TPR_matrix[j,] =medfix_TPR
  medfix_FPR_matrix[j,] =medfix_FPR
  medfix_FDR_matrix[j,] =medfix_FDR
  # bama_TPR_matrix[j,] =bama_TPR
  # bama_FPR_matrix[j,] =bama_FPR
  # bama_FDR_matrix[j,] =bama_FDR
  # mednet_TPR_matrix[j,] =mednet_TPR
  # mednet_FPR_matrix[j,] =mednet_FPR
  # mednet_FDR_matrix[j,] =mednet_FDR
}
colMeans( bama_TPR_matrix)
save(#HIMA_TPR_matrix,
     # HIMA_FPR_matrix,
     # HIMA_FDR_matrix,
     # medfix_TPR_matrix,
     # medfix_FPR_matrix,
     # medfix_FDR_matrix,
     bama_TPR_matrix,
     bama_FPR_matrix,
     bama_FDR_matrix,
     # mednet_TPR_matrix,
     # mednet_FPR_matrix,
     # mednet_FDR_matrix,
     file = '/Users/fuyinghao/Documents/Mediator analysis/Revision/res/bama_nonlinear.RData')







