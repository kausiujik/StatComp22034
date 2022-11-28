#' @title A simulated dataset
#' @name sim_data
#' @description A simulated dataset used to test the individualized regression for clustering subjects
#' @examples
#' \dontrun{
#' data(sim_data)
#' attach(sim_data)
#' result <- indregclu(sim_data, N = 200, K = 2)
#' }
#' @import cluster
#' @import Matrix
#' @importFrom stats rnorm
#' @useDynLib StatComp22034
NULL

#' @title Individualized Regression For Clustering Subjects
#' @description A clustering method based on individualized regression
#' @param data the data frame to be analyzed
#' @param N the number of subjects
#' @param K the number of clustering centers
#' @return a list containing the clustering result and the coefficients for individualized regression
#' @examples
#' \dontrun{
#' data(sim_data)
#' attach(sim_data)
#' result <- indregclu(sim_data, N = 200, K = 2)
#' clu_res <- result$clus_id
#' reg_res <- result$coef
#' }
#' @export
indregclu <- function(data, N, K) {
  update_delta <- function(gamma, beta, nu, lambda, rho) {
    val <- matrix(nrow = K, ncol = N)
    for (k in 1:K) { val[k,] <- (rep(gamma[k], N) - beta - nu / rho)^2 }
    result <- rep(0, N)
    for (i in 1:N) {
      cond <- val[, i] < 2 * lambda / rho
      indT <- which(cond, arr.ind = TRUE)
      if (length(indT) == 0) {result[i] <- beta[i] + nu[i] / rho}
      if (length(indT) == 1) {result[i] <- gamma[indT[1]]}
      if (length(indT) > 1) {result[i] <- gamma[which.min(val[indT, i])]}
    }
    return(result)
  }
  
  y <- data$y; X <- bdiag(diag(data$X)); Z <- data$Z; tildeX <- cbind(X, Z)
  lambda <- 1 / N; rho <- 5
  Q <- Z %*% solve(t(Z) %*% Z) %*% t(Z); P <- diag(N) - Q
  beta0 <- solve(t(X) %*% P %*% X + 0.002 * diag(N)) %*% t(X) %*% P %*% y
  alpha0 <- solve(t(Z) %*% Z) %*% t(Z) %*% (y - X %*% beta0)
  gamma0 <- pam(matrix(beta0, nrow = N), k = K)$medoids
  val <- matrix(nrow = K, ncol = N)
  for (k in 1:K) { val[k, ] <- (beta0@x - gamma0[k,])^2 }
  res_sg <- apply(val, 2, which.min)
  delta0 <- gamma0[res_sg, ] + rnorm(N, 0, 0.1)
  nu0 <- rep(0.1, N)
  
  for (i in 1:1000) {
    betaalpha1 <- solve(t(tildeX) %*% tildeX + bdiag(rho * diag(N), 0)) %*% (t(tildeX) %*% y + c(rho * delta0 - nu0, 0))
    beta1 <- betaalpha1[1:N]
    alpha1 <- betaalpha1[-(1:N)]
    delta1 <- update_delta(gamma0, beta1, nu0, lambda, rho)
    gamma1 <- pam(matrix(delta1, nrow = N), k = K)$medoids
    nu1 <- nu0 + rho * (beta1 - delta1)
    
    res <- beta1 - delta1
    if (sqrt(sum((beta1 - beta0)^2))/N + sqrt(sum((alpha1 - alpha0)^2)) +
        sqrt(sum((gamma1 - gamma0)^2)) < 1e-5 & sqrt(sum(res^2))/N < 1e-4) { break }
    
    beta0 <- beta1
    alpha0 <- alpha1
    delta0 <- delta1
    gamma0 <- gamma1
    nu0 <- nu1
  }
  
  val <- matrix(nrow = K, ncol = N)
  for (k in 1:K) { val[k, ] <- (beta0 - gamma0[k,])^2 }
  res_sg <- apply(val, 2, which.min)
  return(list(clus_id = res_sg, coef = beta0))
}
