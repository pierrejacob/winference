#'@rdname cost_matrix_L1
#'@title cost_matrix_L1
#'@description Compute cost matrix L1 between two matrices of dimension d x N
#'@export
#'
cost_matrix_L1 <- function(x, y){
  return(cost_matrix_L1_(x, y))
}

#'@rdname cost_matrix_L2
#'@title cost_matrix_L2
#'@description Compute cost matrix L2 between two matrices of dimension d x N
#'@export
#'
cost_matrix_L2 <- function(x, y){
  return(cost_matrix_L2_(x, y))
}

#'@rdname cost_matrix_Lp
#'@title cost_matrix_Lp
#'@description Compute cost matrix Lp between two matrices of dimension d x N
#'@export
cost_matrix_Lp <- function(x, y, p){
  return(cost_matrix_Lp_(x, y, p))
}


#'@export
exact_lsap_given_C <- function(w1, w2, C, p = 2){
  solution <- solve_LSAP(C^p)
  return((sum(C[cbind(seq_along(solution), solution)]^p) / nobservations)^(1/p))
}

#'@export
exact_lsap_distance <- function(x1, x2, p = 1, ground_p = 2){
  C <- cost_matrix_Lp(x1, x2, ground_p)
  solution <- solve_LSAP(C^p)
  return((sum(C[cbind(seq_along(solution), solution)]^p) / nobservations)^(1/p))
}

#'@export
get_lsap_to_y <- function(y, p = 1, ground_p = 2){
  library(clue)
  f <- function(z){
    return(exact_lsap_distance(y, z, p, ground_p))
  }
  return(f)
}

#'@export
exact_transport_given_C <- function(w1, w2, C, p = 2){
  a <- transport(w1, w2, costm = C^p, method = "shortsimplex")
  cost <- 0
  for (i in 1:nrow(a)){
    cost <- cost + C[a$from[i], a$to[i]]^p * a$mass[i]
  }
  return(cost^(1/p))
}

#'@export
exact_transport_distance <- function(x1, x2, p = 1, ground_p = 2){
  w1 <- rep(1/ncol(x1), ncol(x1))
  w2 <- rep(1/ncol(x2), ncol(x2))
  C <- cost_matrix_Lp(x1, x2, ground_p)
  # if (ground_p == 1){
  #   C <- cost_matrix_L1(x1, x2)
  # } else {
  #   C <- cost_matrix_L2(x1, x2)
  # }
  library(transport)
  a <- transport(w1, w2, costm = C^p, method = "shortsimplex")
  cost <- 0
  for (i in 1:nrow(a)){
    cost <- cost + C[a$from[i], a$to[i]]^p * a$mass[i]
  }
  return(cost^(1/p))
}

#'@export
get_transport_to_y <- function(y, p = 1, ground_p = 2){
  library(transport)
  f <- function(z){
    return(exact_transport_distance(y, z, p, ground_p))
  }
  return(f)
}

#'@export
sinkhorn_given_C <- function(w1, w2, C, p = 1, eps = 0.05, niterations = 1000){
  epsilon <- eps * median(C^p)
  wass <- winference::wasserstein(w1, w2, C^p, epsilon, niterations)
  ### CORRECTION OF THE MARGINALS
  # explained in the appendix of Coupling of Particle Filters, Jacob Lindsten Schon  (arXiv v2 appendix E)
  Phat <- wass$transportmatrix
  u <- rowSums(Phat)
  utilde <- colSums(Phat)
  alpha <- min(pmin(w1/u, w2/utilde))
  r <- (w1 - alpha * u) / (1 - alpha)
  rtilde <- (w2 - alpha * utilde) / (1 - alpha)
  P <- alpha * Phat + (1 - alpha) * matrix(r, ncol = 1) %*% matrix(rtilde, nrow = 1)
  return(list(uncorrected = (sum(Phat * C^p))^(1/p), corrected = (sum(P * C^p))^(1/p)))
}

#'@export
sinkhorn_distance <- function(x1, x2, p = 1, ground_p = 2, eps = 0.05, niterations = 100){
  w1 <- rep(1/ncol(x1), ncol(x1))
  w2 <- rep(1/ncol(x2), ncol(x2))
  C <- cost_matrix_Lp(x1, x2, ground_p)
  epsilon <- eps * median(C^p)
  wass <- winference::wasserstein(w1, w2, C^p, epsilon, niterations)
  ### CORRECTION OF THE MARGINALS
  # explained in the appendix of Coupling of Particle Filters, Jacob Lindsten Schon  (arXiv v2 appendix E)
  Phat <- wass$transportmatrix
  u <- rowSums(Phat)
  utilde <- colSums(Phat)
  alpha <- min(pmin(w1/u, w2/utilde))
  r <- (w1 - alpha * u) / (1 - alpha)
  rtilde <- (w2 - alpha * utilde) / (1 - alpha)
  P <- alpha * Phat + (1 - alpha) * matrix(r, ncol = 1) %*% matrix(rtilde, nrow = 1)
  return(list(uncorrected = (sum(Phat * C^p))^(1/p), corrected = (sum(P * C^p))^(1/p)))
}

#'@export
get_sinkhorn_to_y <- function(y, p = 1, ground_p = 2, eps = 0.05, niterations = 100){
  f <- function(z){
    return(sinkhorn_distance(y, z, p, ground_p, eps, niterations)$corrected)
  }
  return(f)
}


#'@export
hilbert_distance <- function(x1, x2, p = 1, ground_p = 2, x1sorted = FALSE){
  if (x1sorted){
    orderedx1 <- x1
  } else {
    orderedx1 <- x1[,hilbert_order(x1),drop=F]
  }
  orderedx2 <- x2[,hilbert_order(x2),drop=F]
  distance <- mean(apply(abs(orderedx1 - orderedx2), 2, function(v) (sum(v^ground_p))^(1/ground_p))^p)^(1/p)
  return(distance)
}

#'@export
swap_distance <- function(x1, x2, p = 1, ground_p = 2, init = "hilbert", nsweeps, tolerance){
  if (missing(tolerance) && missing(nsweeps)){
    print("error in swap_distance: needs to specify either tolerance or nsweeps")
    return(Inf)
  }
  if (missing(nsweeps)){
    nsweeps = Inf
  }
  if (missing(tolerance)){
    tolerance <- 0
  }
  Cp <- cost_matrix_Lp(x1, x2, ground_p)^p
  if (init == "hilbert"){
    o1 <- hilbert_order(x1)
    o2 <- hilbert_order(x2)
    permutation <- o2[order(o1)]
  } else {
    permutation <- 1:ncol(x1)
  }
  totalcost <- sum(sapply(1:ncol(x1), function(v) Cp[v,permutation[v]]))
  previous_totalcost <- totalcost
  isweep <- 0
  error <- Inf
  while((isweep < nsweeps)){
    isweep <- isweep + 1
    swapsweep_results <- swapsweep(permutation-1, Cp, totalcost)
    permutation <-  swapsweep_results$permutation + 1
    totalcost <- swapsweep_results$totalcost
    error <- abs(totalcost - previous_totalcost) / ncol(x1)
    if (error < tolerance){
      break
    }
    previous_totalcost <- totalcost
  }
  return(list(distance = (totalcost/ncol(x1))^(1/p), nsweeps = isweep))
}

#'@export
swap_kernel_distance <- function(x1, x2, p = 1, eps = 1, init = "hilbert", nsweeps, tolerance){
  if (missing(tolerance) && missing(nsweeps)){
    print("error in swap_distance: needs to specify either tolerance or nsweeps")
    return(Inf)
  }
  if (missing(nsweeps)){
    nsweeps = Inf
  }
  if (missing(tolerance)){
    tolerance <- 0
  }
  # Cp <- cost_matrix_Lp(x1, x2, ground_p)^p
  Cp <- (1 - exp(-cost_matrix_L2(x1, x2)^2/(2*(eps^2))))^p
  if (init == "hilbert"){
    o1 <- hilbert_order(x1)
    o2 <- hilbert_order(x2)
    permutation <- o2[order(o1)]
  } else {
    permutation <- 1:ncol(x1)
  }
  totalcost <- sum(sapply(1:ncol(x1), function(v) Cp[v,permutation[v]]))
  previous_totalcost <- totalcost
  isweep <- 0
  error <- Inf
  while((isweep < nsweeps)){
    isweep <- isweep + 1
    swapsweep_results <- swapsweep(permutation-1, Cp, totalcost)
    permutation <-  swapsweep_results$permutation + 1
    totalcost <- swapsweep_results$totalcost
    error <- abs(totalcost - previous_totalcost) / ncol(x1)
    if (error < tolerance){
      break
    }
    previous_totalcost <- totalcost
  }
  return(list(distance = (totalcost/ncol(x1))^(1/p), nsweeps = isweep))
}

#'@export
get_hilbert_to_y <- function(y, p = 1, ground_p = 2){
  ysorted <- y[,hilbert_order(y),drop=F]
  f <- function(z){
    return(hilbert_distance(ysorted, z, p, ground_p, x1sorted = TRUE))
  }
  return(f)
}

