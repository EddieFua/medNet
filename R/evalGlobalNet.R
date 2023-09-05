evalGlobalNet <- function(final_net, response, exposure, mediator, test.index, covariate = FALSE, T = 0.6, H = 14, method = "LR") {
  # Evaluate the AUC of the global network generated from `findGlobalNet`
  #
  # Args:
  #   `final_net`: data.frame of the symbolic edges list
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
  #   `auc`: estimated AUC of the network

  # Get number of repeated cross-validation times
  repeat_times <- length(test.index)
  # Get number of folds used in cross-validation
  fold_number <- length(test.index[[1]])

  # Factorize the response for classification instead of regression in RandomForest and SVM
  if (method == "RF" | method == "SVM") {
    response <- as.factor(response)
  }

  list_exposure <- c()
  list_mediator <- c()
  for (i in 1:nrow(final_net)) {
    if(startsWith(final_net[i, 1], "e")) {
      list_exposure <- c(list_exposure, as.numeric(gsub("e", "", final_net[i, 1])))
    } else {
      list_mediator <- c(list_mediator, as.numeric(gsub("m", "", final_net[i, 1])))
    }
    list_mediator <- c(list_mediator, as.numeric(gsub("m", "", final_net[i, 2])))
  }

  list_exposure <- unique(list_exposure)
  list_mediator <- unique(list_mediator)

  temp_comb <- data.frame(response = response, exposure[, list_exposure], mediator[, list_mediator])
  if (class(covariate) != "logical") {
    temp_comb <- cbind(temp_comb, covariate)
  }

  auc <- c()
  flag <- 0
  # Calculate AUC for each repeated set
  for (i in 1:repeat_times) {
    temp_prob <- rep(0, length(response))
    temp_pred <- response
    # Predict for all observations
    for (j in 1:fold_number) {
      test_index0   <- test.index[[i]][[j]]
      testing_set   <- temp_comb[test_index0, ]
      training_set  <- temp_comb[-test_index0, ]
      # Fit classifer
      if (method == "LR") {
        mod <- try(glm(response ~ ., family = binomial, data = training_set, control = list(maxit = 100)))
        if (inherits(mod, 'try-error')) {
          flag <- 1
          break
        }
        temp_prob[test_index0] <- predict(mod, testing_set, type = "response")
        temp_pred[test_index0] <- as.numeric(temp_prob[test_index0] > 0.5)
      } else if (method == "RF") {
        mod <- try(randomForest::randomForest(response ~ ., data = training_set, ntree = 200))
        if (inherits(mod, 'try-error')) {
          flag <- 1
          break
        }
        temp_prob[test_index0] <- predict(mod, testing_set, type = "prob")[, 2]
        temp_pred[test_index0] <- predict(mod, testing_set)
      } else if (method == "SVM") {
        mod <- try(e1071::svm(response ~ ., data = training_set, kernel = "radial", probability = TRUE))
        if (inherits(mod, 'try-error')) {
          flag <- 1
          break
        }
        temp_prob[test_index0] <- attr(predict(mod, testing_set, probability = TRUE), "probabilities")[, 1]
        temp_pred[test_index0] <- predict(mod, testing_set)
      }
    }
    # If failed in fitting model, skip this predictor, and mark it as not predictive
    if (flag == 1) {
      break
    }
    if (as.numeric(pROC::roc(response, temp_prob, direction = "<")$auc)<0.5){
      r = 1-as.numeric(pROC::roc(response, temp_prob, direction = "<")$auc)
      auc <- c(auc, r)
    }else{
      auc <- c(auc,pROC::roc(response, temp_prob, direction = "<")$auc)
    }
  }

  return(mean(auc))
}
