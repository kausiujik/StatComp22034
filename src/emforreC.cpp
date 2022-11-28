# include <Rcpp.h>
using namespace Rcpp;

//' @title EM Algorithm for Repeated Experiment (Rcpp Version)
//' @description Inference on repeated experiment using EM algorithm (Rcpp version).
//' @param data the data frame to be analyzed
//' @return a list containing the estimation results, iteration times and optimal log likelihood.
//' @examples
//' \dontrun{
//' data(reaction_time)
//' attach(reaction_time)
//' result <- emforreC(reaction_time)
//' alpha.hat <- result$alpha; lambda.hat <- result$lambda
//' beta.hat <- result$beta; mu.hat <- result$mu; gamma.hat <- result$gamma
//' sigma.hat <- result$sigma; tau.hat <- result$tau
//' }
//' @export
// [[Rcpp::export]]
List emforreC(DataFrame data) {
  const int n = data.nrows();
  const int t = data.size() - 1;
  NumericVector s = as<NumericVector>(data["type"]);
  NumericMatrix y(n, t);
  for (int j = 1; j <= t; j++) {
    NumericVector yj = as<NumericVector>(data[j]);
    for (int i = 0; i < n; i++) {
      y(i, j-1) = yj[i];
    }
  }
  double sums = 0;
  for (int i = 0; i < n; i++) {
    sums += s[i];
  }
  
  NumericVector alpha0(n), alpha1(n);
  double mu0 = 0, mu1, gamma0 = 100, gamma1;
  double beta0 = 300, beta1, lambda0 = 0.5, lambda1;
  double sigma10 = 10, sigma11, sigma20 = 10, sigma21;
  double tau10 = 10, tau11, tau20 = 10, tau21;
  for (int i = 0; i < n; i++) {
    alpha0[i] = 0;
    for (int j = 0; j < t; j++) {
      alpha0[i] += y(i, j);
    }
    alpha0[i] /= t;
    mu0 += alpha0[i];
  }
  mu0 /= n;
  
  int iters = 0;
  while (TRUE) {
    NumericMatrix xi0(n, t);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < t; j++) {
        double enumij = exp(pow(y(i, j) - alpha0[i] - beta0, 2) / (2 * pow(sigma20, 2)) - pow(y(i, j) - alpha0[i], 2) / (2 * pow(sigma10, 2)));
        xi0(i, j) = s[i] / (1 + enumij * (1 - lambda0) / lambda0 * sigma20 / sigma10);
      }
    }
    
    double sumxi0 = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < t; j++)
        sumxi0 += xi0(i, j);
    lambda1 = sumxi0 / sums / t;
    
    for (int i = 0; i < n; i++) {
      double numer1 = 0, denom1 = 0;
      for (int j = 0; j < t; j++) {
        numer1 += (1 - xi0(i, j)) * y(i, j) / pow(sigma10, 2) + xi0(i, j) * (y(i, j) - beta0) / pow(sigma20, 2);
        denom1 += (1 - xi0(i, j)) / pow(sigma10, 2) + xi0(i, j) / pow(sigma20, 2);
      }
      double numer2 = (1 - s[i]) * mu0 / pow(tau10, 2) + s[i] * (mu0 + gamma0) / pow(tau20, 2);
      double denom2 = (1 - s[i]) / pow(tau10, 2) + s[i] / pow(tau20, 2);
      alpha1[i] = (numer1 + numer2) / (denom1 + denom2);
    }
    
    double wb = 0, ws1 = 0, ws2 = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < t; j++)
        wb += xi0(i, j) * (y(i, j) - alpha1[i]);
    beta1 = wb / sumxi0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < t; j++) {
        ws1 += (1 - xi0(i, j)) * pow(y(i, j) - alpha1[i], 2);
        ws2 += xi0(i, j) * pow(y(i, j) - alpha1[i] - beta1, 2);
      }
      sigma11 = sqrt(ws1 / (n*t - sumxi0));
    sigma21 = sqrt(ws2 / sumxi0);
    
    double wm = 0, wg = 0, wt1 = 0, wt2 = 0;
    for (int i = 0; i < n; i++)
      wm += (1 - s[i]) * alpha1[i];
    mu1 = wm / (n - sums);
    for (int i = 0; i < n; i++)
      wg += s[i] * alpha1[i];
    gamma1 = wg / sums - mu1;
    for (int i = 0; i < n; i++) {
      wt1 += (1 - s[i]) * pow(alpha1[i] - mu1, 2);
      wt2 += s[i] * pow(alpha1[i] - mu1 - gamma1, 2);
    }
    tau11 = sqrt(wt1 / (n - sums));
    tau21 = sqrt(wt2 / sums);
    
    double sumdiffalpha = 0;
    for (int i = 0; i < n; i++)
      sumdiffalpha += abs(alpha0[i] - alpha1[i]);
    
    if (abs(lambda0 - lambda1) < 1e-6 &&
        abs(beta0 - beta1) + abs(mu0 - mu1) + abs(gamma0 - gamma1) < 3e-2 &&
        abs(sigma10 - sigma11) + abs(sigma20 - sigma21) < 2e-4 &&
        abs(tau10 - tau11) + abs(tau20 - tau21) < 2e-4 &&
        sumdiffalpha / n < 1e-2) { break; }
    
    lambda0 = lambda1;
    alpha0 = alpha1;
    beta0 = beta1;
    sigma10 = sigma11;
    sigma20 = sigma21;
    mu0 = mu1;
    gamma0 = gamma1;
    tau10 = tau11;
    tau20 = tau21;
    iters = iters + 1;
  }
  
  NumericVector sigma(2), tau(2);
  sigma[0] = sigma10; sigma[1] = sigma20;
  tau[0] = tau10; tau[1] = tau20;
  
  double like = 0, pi = 3.14159265359;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < t; j++) {
      double dnorm1 = -log(2 * pi * pow(sigma[0], 2)) / 2 - pow(y(i, j) - alpha0[i], 2) / (2 * pow(sigma[0], 2));
      double dnorm2 = -log(2 * pi * pow(sigma[1], 2)) / 2 - pow(y(i, j) - alpha0[i] - beta0, 2) / (2 * pow(sigma[1], 2));
      like += (1 - s[i]) * dnorm1 + s[i] * log((1 - lambda0) * exp(dnorm1) + lambda0 * exp(dnorm2));
    }
    double dnorm1 = -log(2 * pi * pow(tau[0], 2)) / 2 - pow(alpha0[i] - mu0, 2) / (2 * pow(tau[0], 2));
    double dnorm2 = -log(2 * pi * pow(tau[1], 2)) / 2 - pow(alpha0[i] - mu0 - gamma0, 2) / (2 * pow(tau[1], 2));
    like += (1 - s[i]) * dnorm1 + s[i]  * dnorm2;
  }
  
  return List::create(Named("alpha") = alpha0, Named("beta") = beta0, Named("mu") = mu0,
                            Named("gamma") = gamma0, Named("lambda") = lambda0,
                            Named("sigma") = sigma, Named("tau") = tau, Named("iters") = iters,
                            Named("log.lik") = like);
}
