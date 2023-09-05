calAUC = function(k,
                     l,
                     res,
                     response,
                     exposure,
                     mediator,
                     test.index,
                     repeat_times,
                     fold_number,
                     covariate,
                     method){
  temp_comb <- data.frame(response = response, exposure[, k], mediator[, l])
  if (class(covariate) != "logical") {
    temp_comb <- cbind(temp_comb, covariate)
  }

  auc <- c()
  flag <- 0
  # Calculate AUC for each repeated set
  for (i in 1:repeat_times) {
    temp_prob <- rep(NA, length(response))
    temp_pred <- rep(NA, length(response))
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
    auc = as.numeric(pROC::roc(response, temp_prob, direction = "<")$auc)
    res_auc = auc
    res_auc[which(auc<0.5)] = 1-res_auc[which(auc<0.5)]
    res_accuracy = sum(temp_pred == response) / length(response)
    res[i,] = c(res_auc, res_accuracy)
  }
  return(res)
}




