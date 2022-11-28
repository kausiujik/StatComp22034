#' @title A real dataset
#' @name reaction_time
#' @description A real dataset from medical test.
#' @examples
#' \dontrun{
#' data(reaction_time)
#' attach(reaction_time)
#' result <- emforreR(reaction_time)
#' }
#' @useDynLib StatComp22034
NULL

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of R functions (\code{emforreR}) and Rcpp functions (\code{emforreC}).
#' @examples
#' \dontrun{
#' data(reaction_time)
#' attach(reaction_time)
#' tm <- microbenchmark::microbenchmark(
#'   estR = emforreR(reaction_time),
#'   estC = emforreC(reaction_time),
#'   times = 5
#' )
#' print(summary(tm)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm
#' @useDynLib StatComp22034
NULL

#' @title EM Algorithm for Repeated Experiment (R version)
#' @description Inference on repeated experiment using EM algorithm (R version).
#' @param data the data frame to be analyzed
#' @return a list containing the estimation results, iteration times and optimal log likelihood.
#' @examples
#' \dontrun{
#' data(reaction_time)
#' attach(reaction_time)
#' result <- emforreR(reaction_time)
#' alpha.hat <- result$alpha; lambda.hat <- result$lambda
#' beta.hat <- result$beta; mu.hat <- result$mu; gamma.hat <- result$gamma
#' sigma.hat <- result$sigma; tau.hat <- result$tau
#' }
#' @export
emforreR <- function(data) {
  s <- data$type
  n <- length(rownames(data))
  t <- length(data) - 1
  y <- data[colnames(data) != 'type']
  
  log.lik <- function(alpha, beta, mu, gamma, lambda, sigma1, sigma2, tau1, tau2) {
    like <- 0
    for (i in 1:n) {
      for (j in 1:t) {
        like <- like + (1 - s[i]) * log(dnorm(y[i,j], alpha[i], sigma1)) +
          s[i] * log((1 - lambda) * dnorm(y[i,j], alpha[i], sigma1) +
                       lambda * dnorm(y[i,j], alpha[i] + beta, sigma2))
      }
      like <- like + (1 - s[i]) * log(dnorm(alpha[i], mu, tau1)) +
        s[i]  * log(dnorm(alpha[i], mu + gamma, tau2))
    }
    like
  }
  
  alpha0 <- rowMeans(y)
  mu0 <- mean(alpha0)
  gamma0 <- 100; beta0 <- 300; lambda0 <- 0.5
  sigma10 <- sigma20 <- tau10 <- tau20 <- 10
  iters <- 0
  while (TRUE) {
      # update xi
      xi0 <- matrix(nrow = n, ncol = t)
      for (i in 1:n) {
        yi <- as.numeric(y[i,])
        enum <- exp((yi - alpha0[i] - beta0)^2 / (2 * sigma20^2) -
                      (yi - alpha0[i])^2 / (2 * sigma10^2) +
                      log(1 - lambda0) - log(lambda0) +
                      log(sigma20) - log(sigma10))
        xi0[i,] <- s[i] / (1 + enum)
      }
      
      # update lambda
      lambda1 <- sum(xi0) / sum(s) / t
      
      # update alpha
      alpha1 <- numeric(n)
      for (i in 1:n) {
        numer1 <- sum((1 - xi0[i,]) * y[i,] / sigma10^2 + xi0[i,] * (y[i,] - beta0) / sigma20^2)
        denom1 <- sum((1 - xi0[i,]) / sigma10^2 + xi0[i,] / sigma20^2)
        numer2 <- (1 - s[i]) * mu0 / tau10^2 + s[i] * (mu0 + gamma0) / tau20^2
        denom2 <- (1 - s[i]) / tau10^2 + s[i] / tau20^2
        alpha1[i] <- (numer1 + numer2) / (denom1 + denom2)
      }
      
      # update beta
      betai <- numeric(n)
      for (i in 1:n) { betai[i] <- sum(xi0[i,] * (y[i,] - alpha1[i])) }
      beta1 <- sum(betai) / sum(xi0)
      
      # update sigma1
      sigma1i <- numeric(n)
      for (i in 1:n) { sigma1i[i] <- sum((1 - xi0[i,]) * (y[i,] - alpha1[i])^2) }
      sigma11 <- sqrt(sum(sigma1i) / sum(1 - xi0))
      
      # update sigma2
      sigma2i <- numeric(n)
      for (i in 1:n) { sigma2i[i] <- sum(xi0[i,] * (y[i,] - alpha1[i] - beta1)^2) }
      sigma21 <- sqrt(sum(sigma2i) / sum(xi0))
      
      # update mu, gamma, tau1, tau2
      mu1 <- sum((1 - s) * alpha1) / sum(1 - s)
      gamma1 <- sum(s * alpha1) / sum(s) - mu1
      tau11 <- sqrt(sum((1 - s) * (alpha1 - mu1)^2) / sum(1 - s))
      tau21 <- sqrt(sum(s * (alpha1 - mu1 - gamma1)^2) / sum(s))
      
      # break or not?
      if (abs(lambda0 - lambda1) < 1e-6 &&
          mean(abs(beta0 - beta1), abs(mu0 - mu1), abs(gamma0 - gamma1)) < 1e-2 &&
          mean(abs(sigma10 - sigma11), abs(sigma20 - sigma21)) < 1e-4 &&
          mean(abs(tau10 - tau11), abs(tau20 - tau21)) < 1e-4 &&
          mean(abs(alpha0 - alpha1)) < 1e-2) { break }
      
      # iteration
      lambda0 <- lambda1
      alpha0 <- alpha1
      beta0 <- beta1
      sigma10 <- sigma11
      sigma20 <- sigma21
      mu0 <- mu1
      gamma0 <- gamma1
      tau10 <- tau11
      tau20 <- tau21
      iters <- iters + 1
  }
  list(alpha = alpha1, beta = beta1, mu = mu1, gamma = gamma1, lambda = lambda1,
       sigma = c(sigma11, sigma21), tau = c(tau11, tau21), iters = iters,
       log.lik = log.lik(alpha1, beta1, mu1, gamma1, lambda1, sigma11, sigma21, tau11, tau21))
}
