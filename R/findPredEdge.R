library(caret)
library(e1071)
library(energy)
library(pROC)
library(randomForest)
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)


findPredEdge <- function(response, exposure, mediator, test.index, covariate = FALSE, T = 0.6, H = 14, method = "LR", parallization =TRUE, thread = floor(parallel::detectCores()/2)) {
  # Find predictive mediators
  #
  # Args:
  #   `response`: binary numeric vector indicating classifications
  #   `exposure`: data.frame whose rows refer to observations and columns refers to each possible exposure
  #   `mediator`: data.frame whose rows refer to observations and columns refers to each possible mediator
  #   `test.index`: list of index of test observations for each repeat time
  #   `covariate`: data.frame whose rows refer to observations and columns refers to each possible covariate
  #   `T`: threshold used to determine predictivity in each repeat time
  #   `H`: threshold used to determine predictivity in over all repeats
  #   `method`: model used, where `LR` refers to logistic regression
  #                               `RF` refers to random forest
  #                               `SVM` refers to support vector machine
  #
  # Returns:
  #   `result`: a list containing all predictive edges with AUC

  # Get number of repeated cross-validation times
  repeat_times <- length(test.index)
  # Get number of folds used in cross-validation
  fold_number <- length(test.index[[1]])

  # Factorize the response for classification instead of regression in RandomForest and SVM
  if (method == "RF" | method == "SVM") {
    response <- as.factor(response)
  }

  # Calculate AUC for exposures
  result_exposure <- findPredPredictor(response, predictor = exposure, test.index, covariate = covariate, method = method,parallization = parallization)
  # Find indicators for the predictive exposures when there are multiple input exposures
  if (ncol(exposure) > 1) {
    result_exposure$indicator <- (rowSums(result_exposure$auc > T) > H)
  } else {
    result_exposure$indicator <- TRUE
  }
  if (sum(result_exposure$indicator) == 0) {
    cat("No predictive exposure! Try run other exposure sets.\n")
    #cat("No predictive mediator! Try run a single mediator.\n")
    return(0)
    break
  }

  # Calculate the p-values of Brownian distance test between predictive exposure and mediators
  pval <- sapply(which(result_exposure$indicator), function(i) apply(mediator, 2, function(x) energy::dcorT.test(x, exposure[, i])$p.value))
  pval[is.na(pval)] <- 0 #change NA in p-val to 0
  # Find correlated mediators of predictive exposures
  corr_mediator <- mediator[, (rowSums(pval < 0.05) > 0)]
  # Calculate AUC for correlated mediators
  result_mediator <- findPredPredictor(response, corr_mediator, test.index, covariate = covariate, method = method)

  # Create a class of results containing AUC, accuracy, probability and prediction
  result          <- list()
  result$auc      <- array(0, dim = c(ncol(mediator), ncol(exposure), repeat_times))
  result$accuracy <- array(0, dim = c(ncol(mediator), ncol(exposure), repeat_times))

  # Save AUC of exposures for further ranking
  result$exposure_auc <- rowSums(result_exposure$auc) / repeat_times

  res = matrix(NA,nrow = repeat_times,2)
  if (parallization) {
    # Ensure the number of threads specified is not larger than the total cores
    if (parallel::detectCores() < thread) {
      cat("Out of number of cores! Please specify less threads!\n")
    } else {
      mcluster <- parallel::makeCluster(thread, type = "PSOCK")
      doParallel::registerDoParallel(mcluster)
      # Calculate AUC for the pair of each predictive exposure and its correlated mediator using repeated cross validations
      for (k in which(result_exposure$indicator)) {
        iter_num <- which(pval[, which(result_exposure$indicator) == k] < 0.05)
        res_auc <- foreach::foreach(
          l = iter_num,
          .packages = c('medNet','foreach')
        ) %dopar% {
             calAUC(k = k,
                    l = l,
                    res = res,
                    response = response,
                    exposure = exposure,
                    mediator = mediator,
                    test.index,
                    repeat_times,
                    fold_number,
                    covariate,
                    method)
        }
        for (l in which(pval[, which(result_exposure$indicator) == k] < 0.05)) {
          for (i in 1:repeat_times) {
            result$auc[l,k,i] <- res_auc[[which(which(pval[, which(result_exposure$indicator) == k] < 0.05)==l)]][i,1]
            result$accuracy[l,k,i] <- res_auc[[which(which(pval[, which(result_exposure$indicator) == k] < 0.05)==l)]][i,2]}
        }
      }
      stopCluster(mcluster)
    }}else{
      for (k in which(result_exposure$indicator)) {
        for (l in which(pval[, which(result_exposure$indicator) == k] < 0.05)) {
          res_auc = calAUC(k,
                          l,
                          res,
                          response,
                          exposure,
                          mediator,
                          test.index,
                          repeat_times,
                          fold_number,
                          covariate,
                          method)
          for (i in 1:repeat_times) {
            result$auc[l,k,i] <- res_auc[i,1]
            result$accuracy[l,k,i] <- res_auc[i,2]
          }
        }
      }
    }

  # Find predictive mediators for predictive exposure
  result$edge <- matrix(0, length(result_exposure$indicator) * ncol(corr_mediator), 3)
  index <- 1
  for (k in which(result_exposure$indicator)) {
    for (l in which(pval[, which(result_exposure$indicator) == k] < 0.05)) {
      # Check if \sum_r[\hat{AUC}_r(X_i,M_j)>T]>H
      bool1 <- (sum(result$auc[l, k, ] > T) > H)
      # Check if \sum_r[\hat{AUC}_r(X_i,M_j)>\hat{AUC}_r(X_i)]>H
      bool2 <- (sum(result$auc[l, k, ] > result_exposure$auc[k, ]) > H)
      # Check if \sum_r[\hat{AUC}_r(X_i,M_j)>\hat{AUC}_r(M_j)]>H
      bool3 <- (sum(result$auc[l, k, ] > result_mediator$auc[which(which(rowSums(pval < 0.05) > 0) == l), ]) > H)
      # Store the information of predictive edge
      if (bool1 & bool2 & bool3) {
        result$edge[index, 1] <- k
        result$edge[index, 2] <- l
        # result$edge[index, 3] <- sum(result$auc[l, k, ]) / repeat_times #\hat{AUC}_r(X_i,M_j)
        result$edge[index, 3] <- mean(result$auc[l, k, ])
        index <- index + 1
      }
    }
  }
  # Clear redundant rows
  result$edge <- result$edge[(result$edge[, 1] != 0), , drop = FALSE]

  return(result)
}
