## ----eval=F-------------------------------------------------------------------
#  function(data, N, K) {
#    update_delta <- function(gamma, beta, nu, lambda, rho) {
#      val <- matrix(nrow = K, ncol = N)
#      for (k in 1:K) { val[k,] <- (rep(gamma[k], N) - beta - nu / rho)^2 }
#      result <- rep(0, N)
#      for (i in 1:N) {
#        cond <- val[, i] < 2 * lambda / rho
#        indT <- which(cond, arr.ind = TRUE)
#        if (length(indT) == 0) {result[i] <- beta[i] + nu[i] / rho}
#        if (length(indT) == 1) {result[i] <- gamma[indT[1]]}
#        if (length(indT) > 1) {result[i] <- gamma[which.min(val[indT, i])]}
#      }
#      return(result)
#    }
#  
#    y <- data$y; X <- bdiag(diag(data$X)); Z <- data$Z; tildeX <- cbind(X, Z)
#    lambda <- 1 / N; rho <- 5
#    Q <- Z %*% solve(t(Z) %*% Z) %*% t(Z); P <- diag(N) - Q
#    beta0 <- solve(t(X) %*% P %*% X + 0.002 * diag(N)) %*% t(X) %*% P %*% y
#    alpha0 <- solve(t(Z) %*% Z) %*% t(Z) %*% (y - X %*% beta0)
#    gamma0 <- pam(matrix(beta0, nrow = N), k = K)$medoids
#    val <- matrix(nrow = K, ncol = N)
#    for (k in 1:K) { val[k, ] <- (beta0@x - gamma0[k,])^2 }
#    res_sg <- apply(val, 2, which.min)
#    delta0 <- gamma0[res_sg, ] + rnorm(N, 0, 0.1)
#    nu0 <- rep(0.1, N)
#  
#    for (i in 1:1000) {
#      betaalpha1 <- solve(t(tildeX) %*% tildeX + bdiag(rho * diag(N), 0)) %*% (t(tildeX) %*% y + c(rho * delta0 - nu0, 0))
#      beta1 <- betaalpha1[1:N]
#      alpha1 <- betaalpha1[-(1:N)]
#      delta1 <- update_delta(gamma0, beta1, nu0, lambda, rho)
#      gamma1 <- pam(matrix(delta1, nrow = N), k = K)$medoids
#      nu1 <- nu0 + rho * (beta1 - delta1)
#  
#      res <- beta1 - delta1
#      if (sqrt(sum((beta1 - beta0)^2))/N + sqrt(sum((alpha1 - alpha0)^2)) +
#          sqrt(sum((gamma1 - gamma0)^2)) < 1e-5 & sqrt(sum(res^2))/N < 1e-4) { break }
#  
#      beta0 <- beta1
#      alpha0 <- alpha1
#      delta0 <- delta1
#      gamma0 <- gamma1
#      nu0 <- nu1
#    }
#  
#    val <- matrix(nrow = K, ncol = N)
#    for (k in 1:K) { val[k, ] <- (beta0 - gamma0[k,])^2 }
#    res_sg <- apply(val, 2, which.min)
#    return(list(clus_id = res_sg, coef = beta0))
#  }

## ----eval=T, message=F, warning=F---------------------------------------------
library(StatComp22034)
data(sim_data)
attach(sim_data)
result <- indregclu(sim_data, N = 200, K = 2)
result

## ----eval=F-------------------------------------------------------------------
#  function(data) {
#    s <- data$type
#    n <- length(rownames(data))
#    t <- length(data) - 1
#    y <- data[colnames(data) != 'type']
#  
#    log.lik <- function(alpha, beta, mu, gamma, lambda, sigma1, sigma2, tau1, tau2) {
#      like <- 0
#      for (i in 1:n) {
#        for (j in 1:t) {
#          like <- like + (1 - s[i]) * log(dnorm(y[i,j], alpha[i], sigma1)) +
#            s[i] * log((1 - lambda) * dnorm(y[i,j], alpha[i], sigma1) +
#                         lambda * dnorm(y[i,j], alpha[i] + beta, sigma2))
#        }
#        like <- like + (1 - s[i]) * log(dnorm(alpha[i], mu, tau1)) +
#          s[i]  * log(dnorm(alpha[i], mu + gamma, tau2))
#      }
#      like
#    }
#  
#    alpha0 <- rowMeans(y)
#    mu0 <- mean(alpha0)
#    gamma0 <- 100; beta0 <- 300; lambda0 <- 0.5
#    sigma10 <- sigma20 <- tau10 <- tau20 <- 10
#    iters <- 0
#    while (TRUE) {
#        # update xi
#        xi0 <- matrix(nrow = n, ncol = t)
#        for (i in 1:n) {
#          yi <- as.numeric(y[i,])
#          enum <- exp((yi - alpha0[i] - beta0)^2 / (2 * sigma20^2) -
#                        (yi - alpha0[i])^2 / (2 * sigma10^2) +
#                        log(1 - lambda0) - log(lambda0) +
#                        log(sigma20) - log(sigma10))
#          xi0[i,] <- s[i] / (1 + enum)
#        }
#  
#        # update lambda
#        lambda1 <- sum(xi0) / sum(s) / t
#  
#        # update alpha
#        alpha1 <- numeric(n)
#        for (i in 1:n) {
#          numer1 <- sum((1 - xi0[i,]) * y[i,] / sigma10^2 + xi0[i,] * (y[i,] - beta0) / sigma20^2)
#          denom1 <- sum((1 - xi0[i,]) / sigma10^2 + xi0[i,] / sigma20^2)
#          numer2 <- (1 - s[i]) * mu0 / tau10^2 + s[i] * (mu0 + gamma0) / tau20^2
#          denom2 <- (1 - s[i]) / tau10^2 + s[i] / tau20^2
#          alpha1[i] <- (numer1 + numer2) / (denom1 + denom2)
#        }
#  
#        # update beta
#        betai <- numeric(n)
#        for (i in 1:n) { betai[i] <- sum(xi0[i,] * (y[i,] - alpha1[i])) }
#        beta1 <- sum(betai) / sum(xi0)
#  
#        # update sigma1
#        sigma1i <- numeric(n)
#        for (i in 1:n) { sigma1i[i] <- sum((1 - xi0[i,]) * (y[i,] - alpha1[i])^2) }
#        sigma11 <- sqrt(sum(sigma1i) / sum(1 - xi0))
#  
#        # update sigma2
#        sigma2i <- numeric(n)
#        for (i in 1:n) { sigma2i[i] <- sum(xi0[i,] * (y[i,] - alpha1[i] - beta1)^2) }
#        sigma21 <- sqrt(sum(sigma2i) / sum(xi0))
#  
#        # update mu, gamma, tau1, tau2
#        mu1 <- sum((1 - s) * alpha1) / sum(1 - s)
#        gamma1 <- sum(s * alpha1) / sum(s) - mu1
#        tau11 <- sqrt(sum((1 - s) * (alpha1 - mu1)^2) / sum(1 - s))
#        tau21 <- sqrt(sum(s * (alpha1 - mu1 - gamma1)^2) / sum(s))
#  
#        # break or not?
#        if (abs(lambda0 - lambda1) < 1e-6 &&
#            mean(abs(beta0 - beta1), abs(mu0 - mu1), abs(gamma0 - gamma1)) < 1e-2 &&
#            mean(abs(sigma10 - sigma11), abs(sigma20 - sigma21)) < 1e-4 &&
#            mean(abs(tau10 - tau11), abs(tau20 - tau21)) < 1e-4 &&
#            mean(abs(alpha0 - alpha1)) < 1e-2) { break }
#  
#        # iteration
#        lambda0 <- lambda1
#        alpha0 <- alpha1
#        beta0 <- beta1
#        sigma10 <- sigma11
#        sigma20 <- sigma21
#        mu0 <- mu1
#        gamma0 <- gamma1
#        tau10 <- tau11
#        tau20 <- tau21
#        iters <- iters + 1
#    }
#    list(alpha = alpha1, beta = beta1, mu = mu1, gamma = gamma1, lambda = lambda1,
#         sigma = c(sigma11, sigma21), tau = c(tau11, tau21), iters = iters,
#         log.lik = log.lik(alpha1, beta1, mu1, gamma1, lambda1, sigma11, sigma21, tau11, tau21))
#  }

## ----eval=F-------------------------------------------------------------------
#  List emforreC(DataFrame data) {
#    const int n = data.nrows();
#    const int t = data.size() - 1;
#    NumericVector s = as<NumericVector>(data["type"]);
#    NumericMatrix y(n, t);
#    for (int j = 1; j <= t; j++) {
#      NumericVector yj = as<NumericVector>(data[j]);
#      for (int i = 0; i < n; i++) {
#        y(i, j-1) = yj[i];
#      }
#    }
#    double sums = 0;
#    for (int i = 0; i < n; i++) {
#      sums += s[i];
#    }
#  
#    NumericVector alpha0(n), alpha1(n);
#    double mu0 = 0, mu1, gamma0 = 100, gamma1;
#    double beta0 = 300, beta1, lambda0 = 0.5, lambda1;
#    double sigma10 = 10, sigma11, sigma20 = 10, sigma21;
#    double tau10 = 10, tau11, tau20 = 10, tau21;
#    for (int i = 0; i < n; i++) {
#      alpha0[i] = 0;
#      for (int j = 0; j < t; j++) {
#        alpha0[i] += y(i, j);
#      }
#      alpha0[i] /= t;
#      mu0 += alpha0[i];
#    }
#    mu0 /= n;
#  
#    int iters = 0;
#    while (TRUE) {
#      NumericMatrix xi0(n, t);
#      for (int i = 0; i < n; i++) {
#        for (int j = 0; j < t; j++) {
#          double enumij = exp(pow(y(i, j) - alpha0[i] - beta0, 2) / (2 * pow(sigma20, 2)) - pow(y(i, j) - alpha0[i], 2) / (2 * pow(sigma10, 2)));
#          xi0(i, j) = s[i] / (1 + enumij * (1 - lambda0) / lambda0 * sigma20 / sigma10);
#        }
#      }
#  
#      double sumxi0 = 0;
#      for (int i = 0; i < n; i++)
#        for (int j = 0; j < t; j++)
#          sumxi0 += xi0(i, j);
#      lambda1 = sumxi0 / sums / t;
#  
#      for (int i = 0; i < n; i++) {
#        double numer1 = 0, denom1 = 0;
#        for (int j = 0; j < t; j++) {
#          numer1 += (1 - xi0(i, j)) * y(i, j) / pow(sigma10, 2) + xi0(i, j) * (y(i, j) - beta0) / pow(sigma20, 2);
#          denom1 += (1 - xi0(i, j)) / pow(sigma10, 2) + xi0(i, j) / pow(sigma20, 2);
#        }
#        double numer2 = (1 - s[i]) * mu0 / pow(tau10, 2) + s[i] * (mu0 + gamma0) / pow(tau20, 2);
#        double denom2 = (1 - s[i]) / pow(tau10, 2) + s[i] / pow(tau20, 2);
#        alpha1[i] = (numer1 + numer2) / (denom1 + denom2);
#      }
#  
#      double wb = 0, ws1 = 0, ws2 = 0;
#      for (int i = 0; i < n; i++)
#        for (int j = 0; j < t; j++)
#          wb += xi0(i, j) * (y(i, j) - alpha1[i]);
#      beta1 = wb / sumxi0;
#      for (int i = 0; i < n; i++)
#        for (int j = 0; j < t; j++) {
#          ws1 += (1 - xi0(i, j)) * pow(y(i, j) - alpha1[i], 2);
#          ws2 += xi0(i, j) * pow(y(i, j) - alpha1[i] - beta1, 2);
#        }
#        sigma11 = sqrt(ws1 / (n*t - sumxi0));
#      sigma21 = sqrt(ws2 / sumxi0);
#  
#      double wm = 0, wg = 0, wt1 = 0, wt2 = 0;
#      for (int i = 0; i < n; i++)
#        wm += (1 - s[i]) * alpha1[i];
#      mu1 = wm / (n - sums);
#      for (int i = 0; i < n; i++)
#        wg += s[i] * alpha1[i];
#      gamma1 = wg / sums - mu1;
#      for (int i = 0; i < n; i++) {
#        wt1 += (1 - s[i]) * pow(alpha1[i] - mu1, 2);
#        wt2 += s[i] * pow(alpha1[i] - mu1 - gamma1, 2);
#      }
#      tau11 = sqrt(wt1 / (n - sums));
#      tau21 = sqrt(wt2 / sums);
#  
#      double sumdiffalpha = 0;
#      for (int i = 0; i < n; i++)
#        sumdiffalpha += abs(alpha0[i] - alpha1[i]);
#  
#      if (abs(lambda0 - lambda1) < 1e-6 &&
#          abs(beta0 - beta1) + abs(mu0 - mu1) + abs(gamma0 - gamma1) < 3e-2 &&
#          abs(sigma10 - sigma11) + abs(sigma20 - sigma21) < 2e-4 &&
#          abs(tau10 - tau11) + abs(tau20 - tau21) < 2e-4 &&
#          sumdiffalpha / n < 1e-2) { break; }
#  
#      lambda0 = lambda1;
#      alpha0 = alpha1;
#      beta0 = beta1;
#      sigma10 = sigma11;
#      sigma20 = sigma21;
#      mu0 = mu1;
#      gamma0 = gamma1;
#      tau10 = tau11;
#      tau20 = tau21;
#      iters = iters + 1;
#    }
#  
#    NumericVector sigma(2), tau(2);
#    sigma[0] = sigma10; sigma[1] = sigma20;
#    tau[0] = tau10; tau[1] = tau20;
#  
#    double like = 0, pi = 3.14159265359;
#    for (int i = 0; i < n; i++) {
#      for (int j = 0; j < t; j++) {
#        double dnorm1 = -log(2 * pi * pow(sigma[0], 2)) / 2 - pow(y(i, j) - alpha0[i], 2) / (2 * pow(sigma[0], 2));
#        double dnorm2 = -log(2 * pi * pow(sigma[1], 2)) / 2 - pow(y(i, j) - alpha0[i] - beta0, 2) / (2 * pow(sigma[1], 2));
#        like += (1 - s[i]) * dnorm1 + s[i] * log((1 - lambda0) * exp(dnorm1) + lambda0 * exp(dnorm2));
#      }
#      double dnorm1 = -log(2 * pi * pow(tau[0], 2)) / 2 - pow(alpha0[i] - mu0, 2) / (2 * pow(tau[0], 2));
#      double dnorm2 = -log(2 * pi * pow(tau[1], 2)) / 2 - pow(alpha0[i] - mu0 - gamma0, 2) / (2 * pow(tau[1], 2));
#      like += (1 - s[i]) * dnorm1 + s[i]  * dnorm2;
#    }
#  
#    return List::create(Named("alpha") = alpha0, Named("beta") = beta0, Named("mu") = mu0,
#                              Named("gamma") = gamma0, Named("lambda") = lambda0,
#                              Named("sigma") = sigma, Named("tau") = tau, Named("iters") = iters,
#                              Named("log.lik") = like);
#  }

## ----eval=T, message=F, warning=F---------------------------------------------
data(reaction_time)
attach(reaction_time)
resultR <- emforreR(reaction_time)
resultR
resultC <- emforreC(reaction_time)
resultC

## ----eval=T-------------------------------------------------------------------
data(reaction_time)
attach(reaction_time)
tm <- microbenchmark::microbenchmark(
  estR = emforreR(reaction_time),
  estC = emforreC(reaction_time),
  times = 5
)
knitr::kable(summary(tm)[, c(1,3,5,6)])

