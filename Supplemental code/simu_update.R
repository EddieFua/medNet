library(igraph)
StimulateData_update <- function(n_exposure, n_mediator, N, ba_m, ba_power, multi_rate, poisson_rate,
                                 max_mediator_link,true_exposure,num_link){
  # n_feature/n_meta/N: the number of features/meta/sample
  # ba_m: the average number of edges from exposures
  # ba_m/ba_power: the parameters in igraph to simulate the network of metabolites
  # multi_rate: the ratio of features having multiple matching to metabolites
  # poisson_rate: the parameter of Poisson distribution when determing the number of multiple matching
  # max_meta_link/max_feature_link: the cap of number of matchings connected to one metabolite/feature
  # cov_base: the parameters for generating covariance matrix for X
  # true_exposure: the number of true exposures
  # num_link: the number of direct link of single exposure
  ###generate global mediator network
  # n_exposure = 3
  # n_mediator = 500
  # true_exposure = 1
  # N = 400
  #   ba_m = 3
  #   ba_power = 0.5
  #   multi_rate = 0.8
  #   poisson_rate = 1
  #   max_mediator_link = 10
  # true_exposure = true_exposure
  # num_link = 10
  link_num = rep(num_link,true_exposure)
  ba_simulation <- function(num_nodes, ba_m, ba_power){
    graph <- igraph::barabasi.game(n = num_nodes, m = ba_m, power = ba_power, out.pref = F, directed = F, zero.appeal = 0)
    dg <- igraph::degree(graph)
    dist <- igraph::shortest.paths(graph)
    res <- list(graph = graph,dg = dg,dist = dist)
    return(res)
  }
  ba_result <- ba_simulation(n_mediator, ba_m, ba_power)
  exposure_mediator <- matrix(0,n_exposure,n_mediator)
  connected_mediator <- c()  ##save the mediators which are connected with exposures
  functional_link <- matrix(ncol = 2) ##save the link which are between mediators
  mediation_link <- matrix(ncol = 2) ##save the link which are between exposures and mediators
  for(i in 1:true_exposure){
    pool <- 1:n_mediator
    idx1 <- sample(pool, link_num[i],replace = FALSE)
    mediation_link <- rbind(mediation_link,matrix(c(rep(i,length(idx1)),idx1),ncol = 2))
    exposure_mediator[i,idx1] = 1
    for (j in 1:link_num[i]) {
      connected_mediator <- append(connected_mediator,idx1[j])
      idx1_j_neighbor <- which(ba_result$dist[idx1[j],] == (min(ba_result$dist[idx1[j],])+1))   ###find ith exposure's jth mediator's neighbor
      mediator_link <- colSums(exposure_mediator)
      idx1_j_neighbor <- idx1_j_neighbor[which(mediator_link[idx1_j_neighbor] < max_mediator_link)]
      t <- (rpois(1,poisson_rate)+1) * (multi_rate > runif(1))
      if(multi_rate>runif(1)){
        idx2 <- sample(idx1_j_neighbor,min(t,length(idx1_j_neighbor)))
        exposure_mediator[i,idx2] <- 1
        functional_link <- rbind(functional_link,matrix(c(rep(idx1[j],length(idx2)),idx2),ncol = 2))
      }
    }
  }
  mediation_link <- mediation_link[-1,]
  functional_link <- functional_link[-1,]
  link = list(mediation_link = mediation_link,functional_link =functional_link)
  tri_mediator <- Matrix::triu(ba_result$dist)
  mediator_network <- Matrix::which(tri_mediator==1,arr.ind = T)
  ##generate exposure
  res_all = list()
  for (i in 1:3){
    Sigma_matrix_e <- diag(n_exposure)
    diag(Sigma_matrix_e)<-1
    exposure_my <- MASS::mvrnorm(N[[i]], mu = rep(0,n_exposure), Sigma = Sigma_matrix_e)
    Sigma_matrix_m <- 0.3^(ba_result$dist^2)
    Sigma_matrix_m <- corpcor::make.positive.definite(Sigma_matrix_m, tol=1e-10)
    # Sigma_matrix_m <- diag(n_mediator)
    mediator_my <- MASS::mvrnorm(N[[i]], mu = rep(0,n_mediator), Sigma = Sigma_matrix_m)
    id <- which(exposure_mediator==1,arr.ind = T)
    adjust_m <- mediator_my[,id[,2]]+0.5*exposure_my[,id[,1]]
    # for (i in unique(id[,1])){
    #   adjust_e <- exposure_my[,i] + rowSums(0.5*mediator_my[,id[id[,1]==i,][,2]])
    #   exposure_my[,i] = adjust_e
    # }
    mediator_my[,id[,2]] <- adjust_m
    low<-1
    high<-2
    beta_negative_rate <- 0.5
    beta <- runif(true_exposure+1,low,high)
    negative_index <- sample(1:(true_exposure+1),floor((true_exposure+1)*beta_negative_rate))
    beta[negative_index] <- beta[negative_index]*(-1)
    alpha <- rep(0,ncol(mediator_my))
    alpha[id[,2]] <- runif(length(id[,2]),low,high)
    Z <- cbind(matrix(rep(1,N[[i]])),exposure_my[,1:true_exposure])%*%as.matrix(beta,ncol = 1) + mediator_my%*%as.matrix(alpha,ncol = 1)
    response_my <- 1/(1+exp(-Z+median(Z)))
    response_my <- sapply(1:nrow(mediator_my),function(i) if (response_my[i]< 0.5) {response_my[i] = 0 } else{response_my[i] = 1} )
    res_all[[i]] <- list(mediator_network = mediator_network, exposure = exposure_my, response = response_my,  mediator = mediator_my,exposure_mediator = exposure_mediator,link = link)
  }
  
  for (i in 1:3){
    Sigma_matrix_e <- diag(n_exposure)
    diag(Sigma_matrix_e)<-1
    exposure_my <- MASS::mvrnorm(N[[i]], mu = rep(0,n_exposure), Sigma = Sigma_matrix_e)
    # Sigma_matrix_m <- 0.3^(ba_result$dist^2)
    # Sigma_matrix_m <- corpcor::make.positive.definite(Sigma_matrix_m, tol=1e-10)
    Sigma_matrix_m <- diag(n_mediator)
    mediator_my <- MASS::mvrnorm(N[[i]], mu = rep(0,n_mediator), Sigma = Sigma_matrix_m)
    id <- which(exposure_mediator==1,arr.ind = T)
    adjust_m <- mediator_my[,id[,2]]+0.5*exposure_my[,id[,1]]
    # for (i in unique(id[,1])){
    #   adjust_e <- exposure_my[,i] + rowSums(0.5*mediator_my[,id[id[,1]==i,][,2]])
    #   exposure_my[,i] = adjust_e
    # }
    mediator_my[,id[,2]] <- adjust_m
    low<-1
    high<-2
    beta_negative_rate <- 0.5
    beta <- runif(true_exposure+1,low,high)
    negative_index <- sample(1:(true_exposure+1),floor((true_exposure+1)*beta_negative_rate))
    beta[negative_index] <- beta[negative_index]*(-1)
    alpha <- rep(0,ncol(mediator_my))
    alpha[id[,2]] <- runif(length(id[,2]),low,high)
    Z <- cbind(matrix(rep(1,N[[i]])),exposure_my[,1:true_exposure])%*%as.matrix(beta,ncol = 1) + mediator_my%*%as.matrix(alpha,ncol = 1)
    response_my <- 1/(1+exp(-Z+median(Z)))
    response_my <- sapply(1:nrow(mediator_my),function(i) if (response_my[i]< 0.5) {response_my[i] = 0 } else{response_my[i] = 1} )
    res_all[[i+3]] <- list(mediator_network = mediator_network, exposure = exposure_my, response = response_my,  mediator = mediator_my,exposure_mediator = exposure_mediator,link = link)
  }
  
  cl = 10
  ba_result <- ba_simulation(n_mediator, ba_m, ba_power)
  tri_mediator <- Matrix::triu(ba_result$dist)
  mediator_network <- Matrix::which(tri_mediator==1,arr.ind = T)
  ##generate exposure
  Sigma_matrix_e <- diag(n_exposure)
  diag(Sigma_matrix_e)<-0.5
  cl_res = hclust(as.dist(ba_result$dist),method = 'average')
  cl_res = cutree(cl_res, k = cl)
  exposure_mediator <- matrix(0,n_exposure,n_mediator)
  for (i in 1:true_exposure) {
    idx = which(10<table(cl_res)&table(cl_res)<30) ##original 10 & 30
    id = which(cl_res==(idx[i]))
    exposure_mediator[i,id] = 1
  }
  id <- which(exposure_mediator==1,arr.ind = T)
  for (i in 1:3){
    exposure_my <- MASS::mvrnorm(N[i], mu = rep(0,n_exposure), Sigma = Sigma_matrix_e)
    Sigma_matrix_m <- 0.3^(ba_result$dist^2)
    Sigma_matrix_m <- corpcor::make.positive.definite(Sigma_matrix_m, tol=1e-3)
    mediator_my <- MASS::mvrnorm(N[i], mu = rep(0,n_mediator), Sigma = Sigma_matrix_m)
    adjust_m <- mediator_my[,id[,2]]+0.5*exposure_my[,id[,1]]
    mediator_my[,id[,2]] <- adjust_m
    low<-1
    high<-2
    beta_negative_rate <- 0.5
    beta <- runif(true_exposure+1,low,high)
    negative_index <- sample(1:(true_exposure+1),floor((true_exposure+1)*beta_negative_rate))
    beta[negative_index] <- beta[negative_index]*(-1)
    alpha <- rep(0,ncol(mediator_my))
    alpha[id[,2]] <- runif(length(id[,2]),low,high)
    Z <- cbind(matrix(rep(1,N[i])),exposure_my[,1:true_exposure])%*%as.matrix(beta,ncol = 1) + mediator_my%*%as.matrix(alpha,ncol = 1)
    response_my <- 1/(1+exp(-Z+median(Z)))
    response_my <- sapply(1:nrow(mediator_my),function(i) if (response_my[i]< 0.5) {response_my[i] = 0 } else{response_my[i] = 1} )
    link = id
    res_all[[i+6]] <- list(mediator_network = mediator_network, exposure = exposure_my, response = response_my,  mediator = mediator_my,exposure_mediator = exposure_mediator,link = link)
  }
  return(res_all)
}


