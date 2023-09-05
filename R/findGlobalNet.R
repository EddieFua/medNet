library(abind)
library(caret)
library(doParallel)
library(e1071)
library(energy)
library(foreach)
library(igraph)
library(parallel)
library(pROC)
library(randomForest)
findGlobalNet <- function(response, exposure, mediator, mediator.network, covariate = FALSE, repeat.times = 20, fold.number = 5, T = 0.6, H = 14, method = "LR", directed = FALSE, parallization = TRUE,
                          thread = floor(parallel::detectCores()/2)) {
  # Find predictive network for all exposures
  #
  # Args:
  #   `response`: binary numeric vector indicating classifications
  #   `exposure`: data.frame whose rows refer to observations and columns refers to each possible exposure
  #   `mediator`: data.frame whose rows refer to observations and columns refers to each possible mediator
  #   `mediator_network`: data.frame of undirected network chart, where each row refers to a functional link
  #   `repeat.times`: number of repeated cross-validation
  #   `fold.number`: number of fold used in cross-validation
  #   `T`: threshold used to determine predictivity in each repeat time
  #   `H`: threshold used to determine predictivity in over all repeats
  #   `method`: model used, where `LR` refers to logistic regression
  #                               `RF` refers to random forest
  #                               `SVM` refers to support vector machine
  #   `directed`: boolean variable indicating whether the mediator.network is directed or not.
  #   `parallization`: boolean variable indicating whether do parallelization in processing data
  #
  # Returns:
  #  `result`: a list containing the symbolic edges list, undirected network object and the corresponding final AUC

  # Create cross validation folds, and store index of test observations in a list
  test_index <- lapply(1:repeat.times, function(x) caret::createFolds(response, k = fold.number))

  # Use Sys.time as index for ensuring the consistency of the intermidiate results
  id <- as.integer(Sys.time())

  # Do not use parallelization at this moment
  # if (parallel) {
  #   return("Parallelization under testing and do not use it!")
  # }
  # Do parallelization
  # Find all predictive exposures and its predictive mediators
  result_edge <- findPredEdge(response, exposure, mediator, test_index, covariate = covariate, T = T, H = H, method = method,parallization,thread = thread)
  if (class(result_edge)!="list") {
    cat("No predictive exposure! Try run other exposure sets.\n")
    return(0)
    break
  }

  # Regenerate the predictive exposures
  result_exposure <- findPredPredictor(response, exposure, test_index, covariate = covariate, parallization,method = method,thread = thread)
  result_exposure$indicator <- (rowSums(result_exposure$auc > T) > H)
  result_exposure$rank <- order(rowSums(result_exposure$auc) / repeat.times, decreasing = T)

  # Save the intermidiate results
  save(result_edge, file = paste0("./predictive_edge_", id, ".RData"))
  save(result_exposure, file = paste0("./predictive_exposure_", id, ".RData"))

  # If no exposure found, stop the process
  if(any(result_exposure$indicator)==FALSE) stop("No potential exposure is identified! Try to set a smaller T or H!")

  # Find predictive network for each single predictive exposure
  single_network <- list()
  for (k in unique(result_edge$edge[, 1])) {
    # Substract all the information of a single predictive exposure with its predictive mediators
    single_exposure          <- list()
    single_exposure$exposure <- exposure[, k, drop = FALSE]
    single_exposure$auc      <- result_edge$auc[, k, ]
    single_exposure$accuracy <- result_edge$accuracy[, k, ]
    single_exposure$edge     <- result_edge$edge[result_edge$edge[, 1] == k, , drop = FALSE]
    # If the mediator network is directed
    if (directed) {
      # Find predictive network for single exposure
      single_network[[paste0("e",k)]] <- findSingleNet2(response, single_exposure, mediator, mediator.network, test_index, covariate = covariate)
    } else {
      single_network[[paste0("e",k)]] <- findSingleNet(response, single_exposure, mediator, mediator.network, test_index, covariate = covariate)
      # single_network[[paste0("e",k)]] <- findSingleNet(response, single_exposure, mediator, mediator_network, test_index, covariate = covariate,parallization = parallization)
    }
  }
  # Save the intermidiate results
  save(single_network, file = paste0("./single_network_", id, ".RData"))

  # Load the predictive exposures
  if (file_test("-f", paste0("./predictive_exposure_", id, ".RData"))) {
    load(paste0("./predictive_exposure_", id, ".RData"))
    # Rank the exposure by their estimated AUC
    result_exposure$rank <- order(rowSums(result_exposure$auc) / 20, decreasing = TRUE)
  } else {
    cat("Have you stoped the main function? No intermidiate files identified!\n")
    # In case of inconsistency of the intermidiate files
    result_exposure <- findPredPredictor(response, exposure, test_index, covariate = covariate, method = method)
    result_exposure$indicator <- (rowSums(result_exposure$auc > T) > H)
    result_exposure$rank <- order(rowSums(result_exposure$auc) / 20, decreasing = TRUE)
  }

  # Generate the global network
  names  <- colnames(exposure)
  mnames <- colnames(mediator)

  final_net <- c()
  for (k in unique(result_edge$edge[, 1])) {
    # Attach the mediation links
    final_net <- rbind(final_net, data.frame(from = paste0("e", k), to = paste0("m", single_network[[paste0("e", k)]]$mediation_link[, 2]), type = "mediation_link"))
    # Attach the functional links
    tlink <- single_network[[paste0("e", k)]]$functional_link
    # In case of returning a vector instead of a matrix
    if (class(tlink)[1] == "numeric") {
      tlink = matrix(tlink, ncol = 3)
    }
    # Check whether there is any functional link
    if (nrow(tlink) != 0) {
      final_net <- rbind(final_net, data.frame(from = paste0("m", tlink[, 1]), to = paste0("m", tlink[, 2]), type = "functional_link"))
    }
  }

  final_net[, 1] <- as.character(final_net[, 1])
  final_net[, 2] <- as.character(final_net[, 2])
  # Remove repeated functional links
  final_net <- unique(final_net[, 1:3])

  # Estimate the AUC of the final network
  auc <- evalGlobalNet(final_net, response, exposure, mediator, test_index, covariate = covariate, T = T, H = H, method = method)

  # Create vertice types
  node <- unique(c(final_net[, 1], final_net[, 2]))
  vtype <- rep(2, length(node))
  vtype[grep("^e", node, perl = T)] <- 1

  # Assign the original names of exposures and mediators in the final_net
  for (i in 1:nrow(final_net)) {
    # Use pseudo names if no names assigned for mediators
    if (length(mnames) != 0) {
      if (length(names) != 0 & startsWith(final_net[i, 1], "e")) {
        final_net[i, 1] <- names[as.numeric(gsub("e", "", final_net[i, 1]))]
      } else {
        final_net[i, 1] <- mnames[as.numeric(gsub("m", "", final_net[i, 1]))]
      }
      final_net[i, 2] <- mnames[as.numeric(gsub("m", "", final_net[i, 2]))]
    } else {
      if (length(names) != 0 & startsWith(final_net[i, 1], "e")) {
        final_net[i, 1] <- names[as.numeric(gsub("e", "", final_net[i, 1]))]
      }
    }}
  # Create an undirected graph
  node <- unique(c(final_net[, 1], final_net[, 2]))
  net <- igraph::graph_from_data_frame(final_net, directed = F, vertices = data.frame(node, vtype))

  # Save results
  result           <- list()
  result$AUC       <- auc
  result$net       <- net
  result$final_net <- final_net
  # Visualize the global network
  result$plot      <- plotMedNet(result)

  return(result)
}
