library(caret)
library(e1071)
library(pROC)
library(randomForest)
library(parallel)
library(foreach)
library(doParallel)

findPredPredictor <-function(response,
                                predictor,
                                test.index,
                                covariate = FALSE,
                                method = "LR",
                                classifier = FALSE,
                                parallization = T,
                                thread = floor(parallel::detectCores()/2)) {
  # Find predictive predictors, which can be exposures or mediators
  #
  # Args:
  #   `response`: binary numeric vector indicating classifications
  #   `exposure`: data.frame whose rows refer to observations and columns refers to each possible exposure
  #   `predictor`: data.frame whose rows refer to observations and columns refers to each possible predictor
  #   `test.index`: list of index of test observations for each repeat time
  #   `covariate`: data.frame whose rows refer to observations and columns refers to each possible covariate
  #   `method`: model used, where `LR` refers to logistic regression
  #                               `RF` refers to random forest
  #                               `SVM` refers to support vector machine
  #                               `userdefined` refers to any classifier user want to use
  #   `classifier`: a classifier user want to use when `method` is given as `userdefined`
  #   `thread`:number of cores used in parallel processing if allowed, default is floor(parallel::detectCores()/2)
  # Returns:
  #   `result`: a list of AUC, accuracy, probability and prediction for single predictor

  # Get number of repeated cross-validation times
  repeat_times <- length(test.index)
  # Get number of fold used in cross-validation
  fold_number <- length(test.index[[1]])
  # the function that we propose to parralize

  # Factorize the response for classification instead of regression in RandomForest and SVM
  if (method == "RF" | method == "SVM") {
    response <- as.factor(response)
  }

  # Create a class of results containing AUC, accuracy, probability and prediction
  result          <- list()
  result$auc      <- matrix(0, ncol(predictor), repeat_times)
  result$accuracy <- matrix(0, ncol(predictor), repeat_times)
  my_res = list()
  res = matrix(NA,repeat_times,2)
  l = 0
  if (parallization) {
    # Ensure the number of threads specified is not larger than the total cores
    if (parallel::detectCores() < thread) {
      cat("Out of number of cores! Please specify less threads!\n")
    } else {
      # Register a cluster which creates brand new R Sessions and inherit nothing from the master
      mcluster <- parallel::makeCluster(thread, type = "PSOCK")
      doParallel::registerDoParallel(mcluster)
      n = ncol(predictor)
      my_res = foreach::foreach(
        k = 1:n,
        .packages = c('medNet','foreach')
        # .combine = rbind
      ) %dopar% {
        calAUC(
          k,
          l,
          res,
          response,
          predictor,
          predictor,
          test.index,
          repeat_times,
          fold_number,
          covariate,
          method
        )
      }
      # Stop cluster
      stopCluster(mcluster)
    }}else{
      for (k in 1:ncol(predictor)) {
        my_res[[k]] = calAUC(
          k,
          l,
          res,
          response,
          predictor,
          predictor,
          test.index,
          repeat_times,
          fold_number,
          covariate,
          method
        )
      }}
  for (i in 1:length(my_res)){
    result$auc[i,] = my_res[[i]][,1]
    result$accuracy[i,] = my_res[[i]][,2]
  }
  return(result)
}







