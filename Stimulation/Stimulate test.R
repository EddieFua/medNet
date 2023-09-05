library(igraph)
n_exposure = 1
n_mediator = 1000
N = 700
ba_m = 3
ba_power = 1
multi_rate = 0.5
poisson_rate = 1
max_mediator_link = 10
link_num = 20
ba_simulation <- function(num_nodes, ba_m, ba_power) {
  graph <-
    igraph::barabasi.game(
      n = num_nodes,
      m = ba_m,
      power = ba_power,
      out.pref = F,
      directed = F,
      zero.appeal = 0
    )
  dg <- igraph::degree(graph)
  dist <- igraph::shortest.paths(graph)
  res <- list(graph = graph, dg = dg, dist = dist)
  return(res)
}
ba_result <- ba_simulation(n_mediator, ba_m, ba_power)
exposure_mediator <- matrix(0, n_exposure, n_mediator)
connected_mediator <- c()  ##save the mediators which are connected with exposures
functional_link <- matrix(ncol = 2) ##save the link which are between mediators
for (i in 1:n_exposure) {
  pool <- 1:n_mediator
  idx1 <- sample(pool, link_num[i], replace = FALSE)
  exposure_mediator[i, idx1] = 1
  for (j in 1:link_num[i]) {
    connected_mediator <- append(connected_mediator, idx1[j])
    idx1_j_neighbor <-
      which(ba_result$dist[idx1[j], ] == (min(ba_result$dist[idx1[j], ]) + 1))   ###find ith exposure's jth mediator's neighbor
    mediator_link <- colSums(exposure_mediator)
    idx1_j_neighbor <-
      idx1_j_neighbor[which(mediator_link[idx1_j_neighbor] < max_mediator_link)]
    t <- (rpois(1, poisson_rate) + 1) * (multi_rate > runif(1))
    if (multi_rate > runif(1)) {
      idx2 <- sample(idx1_j_neighbor, min(t, length(idx1_j_neighbor)))
      exposure_mediator[i, idx2] <- 1
      functional_link <-
        rbind(functional_link, matrix(c(rep(
          idx1[j], length(idx2)
        ), idx2), ncol = 2))
    }
  }
}
functional_link <- functional_link[-1, ]
tri_mediator <- Matrix::triu(ba_result$dist)
mediator_network <- Matrix::which(tri_mediator == 1, arr.ind = T)
##generate exposure
Sigma_matrix_e <- diag(n_exposure) * 0.5
exposure_my <- MASS::mvrnorm(N, mu = rep(0, n_exposure), Sigma = Sigma_matrix_e)
Sigma_matrix_m <- diag(n_mediator) * 0.5
mediator_my <- MASS::mvrnorm(N, mu = rep(0, n_mediator), Sigma = Sigma_matrix_m)


single_low = seq(0.1,0.9,0.1)
# single_high = 0.6
low_range = seq(0.1,0.9,0.1)
# high_range = 0.6
coe_range = seq(0.1,1,0.1)
comb_range = matrix(NA,ncol = 5)
for (i in 1:length(single_low)) {
  single_high = ((1:((1-single_low[i])/0.1))-1)*0.1+single_low[i]+0.1
  low_range = ((1:((1-single_low[i])/0.1))-1)*0.1+single_low[i]
  for (j in 1:length(low_range)) {
    high_range = ((1:((1-low_range[j])/0.1))-1)*0.1+low_range[j]+0.1
    comb_range = rbind(comb_range,as.matrix(crossing(low_range[j],high_range,single_low[i],single_high,coe_range)))
  }
}
comb_range = comb_range[-1,]
num = rep(NA,nrow(comb_range))
for (j in 1:nrow(comb_range)) {
  mediator_add <- matrix(0, nrow(mediator_my), ncol(mediator_my))
  for (i in 1:nrow(mediator_network)) {
    coe <- runif(1, comb_range[j,3], comb_range[j,4])
    mediator_add[, mediator_network[i, 2]] = coe * mediator_my[, mediator_network[i, 1]] + mediator_add[, mediator_network[i, 2]]
    mediator_add[, mediator_network[i, 1]] = coe * mediator_my[, mediator_network[i, 2]] + mediator_add[, mediator_network[i, 1]]
  }

  mediator <- mediator_my + mediator_add
  id <- which(exposure_mediator == 1, arr.ind = T)
  coe = comb_range[j,5]
  mediator[, id[, 2]] <- mediator[, id[, 2]] + coe*exposure_my[, id[, 1]]
  low <- comb_range[j,1]
  high <- comb_range[j,2]
  beta_negative_rate <- 0.5
  beta <- runif(n_exposure + 1, low, high)
  negative_index <- sample(1:(n_exposure + 1), floor((n_exposure + 1) * beta_negative_rate))
  beta[negative_index] <- beta[negative_index] * (-1)
  alpha <- rep(0, ncol(mediator))
  alpha[id[, 2]] <- runif(length(id[, 2]), low, high)
  Z <- cbind(matrix(rep(1, N)), exposure_my) %*% beta + mediator %*% alpha
  response_my <- 1 / (1 + exp(-Z))
  # plot(Z,1 / (1 + exp(-Z)), main = paste0('low=',comb_range[i,1],' high=',comb_range[i,2],' single_low=',comb_range[i,3],' single_high=',
  # comb_range[i,4]))
  num[j] = sum(response_my <= 0.9&response_my>=0.1)
}

which(num==max(num))

pdf('/Users/fuyinghao/Documents/medNetmy/res/logit.pdf')
plot(Z,response_my)
dev.off()
