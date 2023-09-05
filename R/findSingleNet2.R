library(caret)
library(e1071)
library(energy)
library(pROC)
library(randomForest)

findSingleNet2 <- function(response, single.exposure, mediator, mediator.network, test.index, covariate = FALSE, T = 0.6, H = 14, method = "LR") {
  # Find predictive network for single exposure with directed mediator network
  #
  # Args:
  #   `response`: binary numeric vector indicating classifications
  #   `single.exposure`: list of information of a single predictive exposure with its predictive mediators
  #   `mediator`: data.frame whose rows refer to observations and columns refers to each possible mediator
  #   `mediator.network`: data.frame of undirected network chart, where each row refers to a functional link
  #   `test.index`: list of index of test observations for each repeat time
  #   `covariate`: data.frame whose rows refer to observations and columns refers to each possible covariate
  #   `T`: threshold used to determine predictivity in each repeat time
  #   `H`: threshold used to determine predictivity in over all repeats
  #   `method`: model used, where `LR` refers to logistic regression
  #                               `RF` refers to random forest
  #                               `SVM` refers to support vector machine
  #
  # Returns:
  #   `result`: a list of mediation links and functional links of a single exposure

  # Get number of repeated cross-validation times
  repeat_times <- length(test.index)
  # Get number of fold used in cross-validation
  fold_number <- length(test.index[[1]])
  # Get the value of exposure
  exposure <- single.exposure$exposure

  # Factorize the response for classification instead of regression in RandomForest and SVM
  if (method == "RF" | method == "SVM") {
    response <- as.factor(response)
  }

  # Start from the predictive mediator with largest AUC
  edge <- single.exposure$edge[order(single.exposure$edge[, 3], decreasing = TRUE), , drop = FALSE]

  # Record whether a functional link in mediator network has been visited or not
  visit_record <- rep(FALSE, nrow(mediator.network))

  # Save found C=\{X_i,M_{j_1},...,M_{j_n}\}
  final_edge <- matrix(0, nrow(mediator.network), 3)
  index <- 1

  for (l in edge[, 2]) {
    comb <- data.frame(response = response, exposure = exposure, mediator[, l])
    colnames(comb)[ncol(comb)] <- paste0("mediator_", l)
    if (class(covariate) != "logical") {
      comb <- cbind(comb, covariate)
    }

    # Find the functional links in the mediator network involving the predictive mediator_l
    functional_link <- which((rowSums(mediator.network == l) == 1) & !visit_record)
    # Save the nodes we found
    list <- c(l)

    # Initiate \hatAUC(C)
    current_auc <- single.exposure$auc[l, ]

    # Go through all possible functional links
    while (length(functional_link) > 0) {
      temp_auc <- matrix(0, length(functional_link), repeat_times)
      for (r in functional_link) {
        # Find mediator_m linked to the predictive mediator mediator_l
        m <- mediator.network[r, !(mediator.network[r, ] %in% list)]
        temp_comb <- cbind(comb, mediator[, m])
        colnames(temp_comb)[ncol(temp_comb)] <- paste0("mediator_", m)
        # Annotate whether there is error in fitting classifier
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
          temp_auc[which(functional_link == r), i] <- pROC::roc(response, temp_prob, direction = "<")$auc
          # result$accuracy[k, i] <- sum(temp_pred == response) / length(response)
        }
      }
      # Check if \sum_r[\hat{AUC}_r(C,Ne_i)>\hat{AUC}_r(C)]>H
      bool <- rowSums(temp_auc > current_auc) > H
      if (sum(bool) > 0) {
        # Save current AUC in case of not added
        temp_current_auc <- current_auc
        # Calculate AUC for the significant Ne_i
        current_auc <- temp_auc
        temp_auc <- rowSums(temp_auc) / repeat_times

        # Save argmax_{Ne_i}
        argmax <- which(temp_auc == max(temp_auc[bool]))
        # Indicate whether any directed functional link is found
        added <- FALSE

        # mediator_argmax may appear several times since it connects to several predictive mediators in the `list`
        for (argmax_i in argmax) {
          # Only functional links from col1 to col2 will be recorded
          if (mediator.network[functional_link[argmax_i], 1] == l) {
            final_edge[index, 1] <- mediator.network[functional_link[argmax_i], 1]
            final_edge[index, 2] <- mediator.network[functional_link[argmax_i], 2]
            final_edge[index, 3] <- max(temp_auc[bool])
            index <- index + 1
            current_auc <- current_auc[argmax_i, ]
            added <- TRUE
          }
        }

        # Mark redundant functional link and argmax functional link as visited
        visit_record[functional_link[!bool]]  <- TRUE
        visit_record[functional_link[argmax]] <- TRUE

        # Remove redundant functional link and argmax functional link
        bool[argmax] <- FALSE
        argmax <- mediator.network[functional_link[argmax[1]], ]
        argmax <- as.numeric(argmax[!(argmax %in% list)])
        functional_link <- functional_link[bool]

        # If found any directed functional link, keep going
        if (added) {
          # Add the functional links in the mediator network involving the predictive mediator_argmax
          functional_link <- c(functional_link, which((rowSums(mediator.network == argmax) == 1) & !visit_record))
          list <- c(list, argmax)
          # Move to find the next
          comb <- cbind(comb, mediator[, argmax])
          colnames(comb)[ncol(comb)] <- paste0("mediator_", argmax)
        } else {
          # Restore the previous AUC for not added case
          current_auc <- temp_current_auc
        }
      } else {
        # Mark redundant functional link as visited
        visit_record[functional_link[!bool]]  <- TRUE
        # Remove redundant functional link
        functional_link <- functional_link[bool]
      }
    }
  }

  # Save result
  final_edge <- final_edge[final_edge[, 1] != 0, , drop = FALSE]
  result <- list()
  result$mediation_link  <- edge
  result$functional_link <- final_edge

  return(result)
}
