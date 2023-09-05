library(igraph)
StimulateData <- function(n_exposure, n_mediator, N,link_num, ba_m, ba_power, multi_rate, poisson_rate,
                          max_mediator_link){
  # n_feature/n_meta/N: the number of features/meta/sample
  # ba_m: the average number of edges from exposures
  # ba_m/ba_power: the parameters in igraph to simulate the network of metabolites
  # multi_rate: the ratio of features having multiple matching to metabolites
  # poisson_rate: the parameter of Poisson distribution when determing the number of multiple matching
  # max_meta_link/max_feature_link: the cap of number of matchings connected to one metabolite/feature
  # cov_base: the parameters for generating covariance matrix for X
  # n_exposure = 1
  # n_mediator = 5000
  # N = 1000
  # ba_m = 3
  # ba_power = 1
  # multi_rate = 0.5
  # poisson_rate = 0.5
  # max_mediator_link = 5
  # link_num = 20
  ###generate global mediator network
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
  for(i in 1:n_exposure){
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
  Sigma_matrix_e <- diag(n_exposure)*0.5
  exposure_my <- MASS::mvrnorm(N, mu = rep(0,n_exposure), Sigma = Sigma_matrix_e)
  Sigma_matrix_m <- diag(n_mediator)*0.5
  mediator_my <- MASS::mvrnorm(N, mu = rep(0,n_mediator), Sigma = Sigma_matrix_m)
  mediator_add <- matrix(0,nrow(mediator_my),ncol(mediator_my))
  for (i in 1:nrow(mediator_network)) {
    # coe <- runif(1,0.35,0.6)
    coe <- runif(1,0.1,0.3)
    mediator_add[, mediator_network[i, 2]] = coe * mediator_my[, mediator_network[i, 1]] + mediator_add[, mediator_network[i, 2]]
    mediator_add[, mediator_network[i, 1]] = coe * mediator_my[, mediator_network[i, 2]] + mediator_add[, mediator_network[i, 1]]
  }
  mediator_my <- mediator_my+mediator_add
  id <- which(exposure_mediator==1,arr.ind = T)
  mediator_my[,id[,2]] <- mediator_my[,id[,2]]+0.1*exposure_my[,id[,1]]

  # low <- 0.25
  low<-0.1
  high<-0.2
  # high <- 0.6
  beta_negative_rate <- 0.5
  beta <- runif(n_exposure+1,low,high)
  negative_index <- sample(1:(n_exposure+1),floor((n_exposure+1)*beta_negative_rate))
  beta[negative_index] <- beta[negative_index]*(-1)
  alpha <- rep(0,ncol(mediator_my))
  alpha[id[,2]] <- runif(length(id[,2]),low,high)
  Z <- cbind(matrix(rep(1,N)),exposure_my)%*%beta + mediator_my%*%alpha
  response_my <- 1/(1+exp(-Z))
  response_my <- sapply(1:nrow(mediator_my),function(i) if (response_my[i]< 0.5) {response_my[i] = 0 } else{response_my[i] = 1} )
  res <- list(mediator_network = mediator_network, exposure = exposure_my, response = response_my,  mediator = mediator_my,exposure_mediator = exposure_mediator,link = link)
  return(res)
}


