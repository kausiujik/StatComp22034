## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----fig.height=4, fig.width=6, comment=""------------------------------------
attach(iris)  # "iris" is a famous dataset
colnames(iris)
plot(Sepal.Length, Petal.Length)

## ----fig.height=4, fig.width=6------------------------------------------------
scatterplot3d::scatterplot3d(Sepal.Length, Sepal.Width, Petal.Length, highlight.3d = TRUE)

## ----fig.height=4, fig.width=6, message=FALSE---------------------------------
library(ggplot2)
library(grid)
ggplot(iris, aes(Sepal.Length, Petal.Length)) +
         geom_smooth(method = "lm") +  # Use linear regression to fit the model.
         geom_point() +
         facet_wrap(~Species)  # Under different species.

## ---- comment = ""------------------------------------------------------------
knitr::kable(head(iris), align = "l")
xtable::xtable(head(iris))

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----fig.height=4, fig.width=6------------------------------------------------
set.seed(22034)
pareto <- function(a,b) {
  # Generate random numbers with cdf F(x)
  u <- runif(10000)
  x <- b*(1-u)^(-1/a)
  
  # Draw the histogram of random numbers generated
  hist(x, breaks = "scott", prob = TRUE, main = paste('Pareto(',a,',',b,')'))
  
  # Draw the density function f(x)
  y <- seq(0, max(x), 0.1)
  lines(y, a*b^a/(y^(a+1)))
}

pareto(2, 2)

## ----fig.height=4, fig.width=6------------------------------------------------
set.seed(22034)
beta <- function(a,b) {
  # Calculate constant c
  x0 <- (a-1)/(a+b-2)
  c <- x0^(a-1)*(1-x0)^(b-1)  # constant in pdf can be ignored
  
  # Generate random numbers with pdf f(x)
  n <- 1000
  k <- 0
  y <- numeric(n)
  while (k < n) {
    u <- runif(1)
    x <- runif(1) # random variate from g(x)
    if (x^(a-1)*(1-x)^(b-1) / c > u) {
      # accept x
      k <- k + 1
      y[k] <- x
    }
  }
  
  # Draw the histogram of random numbers generated
  hist(y, breaks = "scott", prob = TRUE, main = paste('Beta(',a,',',b,')'), xlab = "x")
  
  # Draw the density function f(x)
  z <- seq(0, 1, 0.01)
  lines(z, z^(a-1)*(1-z)^(b-1)*gamma(a+b)/gamma(a)/gamma(b))
}

beta(3, 2)

## -----------------------------------------------------------------------------
set.seed(22034)
expgamma <- function(r, beta) {
  # Generate random numbers from the mixture
  n <- 1000
  x <- rgamma(n, r, beta)
  y <- rexp(n, x)
  return(y)
}

r <- 4; beta <- 2
rnd <- expgamma(r, beta)

## ----fig.height=4, fig.width=6------------------------------------------------
# Draw the histogram of random numbers generated
hist(rnd, breaks = "scott", prob = TRUE, main = paste('Pareto(',r,',',beta,')'), xlab = "y")

# Draw the density function f(y)
y <- seq(0, max(rnd), 0.01)
lines(y, r*beta^r/(beta+y)^(r+1))

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----fig.height=4, fig.width=6------------------------------------------------
quick_sort <- function(x) {
  num <- length(x)
  if (num==0 || num==1) {return(x)}
  else {
    a <- x[1]
    y <- x[-1]
    lower <- y[y<a]
    upper <- y[y>=a]
    return (c(quick_sort(lower), a, quick_sort(upper)))
  }
}

n <- c(1e4, 2e4, 4e4, 6e4, 8e4)
set.seed(22034)
an <- numeric(5)
for (i in 1:5) {
  numbers <- sample(1:n[i])
  an[i] <- mean(replicate(100, system.time(quick_sort(numbers))[1]))
}

tn <- n*log(n)
lm.fit <- lm(an ~ tn)
plot(tn, an, col='blue', pch=2)
abline(lm.fit, col='red', lwd=2)

## -----------------------------------------------------------------------------
mc <- function(m) {
  x <- runif(m)
  return(mean(exp(x)))
}

## -----------------------------------------------------------------------------
anti <- function(m) {
  x <- runif(m/2)
  y <- 1 - x
  return((mean(exp(x)) + mean(exp(y)))/2)
}

## ----comment=''---------------------------------------------------------------
iters <- 1000
m <- 10000
set.seed(22034)
theta1 <- theta2 <- numeric(iters)
for (i in 1:iters) {  # Iteration for times
  theta1[i] <- mc(m)
  theta2[i] <- anti(m)
}
thetahat1 <- mean(theta1)
thetahat2 <- mean(theta2)
c(thetahat1, thetahat2)
varhat1 <- var(theta1)
varhat2 <- var(theta2)
c(varhat1, varhat2)
(varhat1 - varhat2) / varhat1

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----fig.height=4, fig.width=6------------------------------------------------
g <- function(x) {x^2 * exp(-x^2 / 2) / sqrt(2 * pi) * (x>=1)}
plot(g, xlim = c(1, 5))

## ----fig.height=6, fig.width=6------------------------------------------------
f1 <- function(x) {exp(1 - x) * (x>=1)}
f2 <- function(x) {2 * exp(-(x - 1)^2 / 2) / sqrt(2 * pi) * (x>=1)}
x <- seq(1, 5, 0.01)
plot(g, xlim = c(1, 5), ylim = c(0, 1), lwd = 2)
lines(x, f1(x), lty = 2, col = 2, lwd = 2)
lines(x, f2(x), lty = 3, col = 3, lwd = 2)
legend("topright", c(expression(g(x)==e^{-x^2/2}/sqrt(2*pi)),
                     expression(f1(x)==e^{-x+1}),
                     expression(f2(x)==2*e^{-(x-1)^2/2}/sqrt(2*pi))),
       inset = 0.02, lty = 1:3, col = 1:3, lwd = 2)

## ----comment=''---------------------------------------------------------------
c(integrate(f1, 1, Inf)$value, integrate(f2, 1, Inf)$value)

## ----fig.height=6, fig.width=6------------------------------------------------
x <- seq(1, 5, 0.01)
gf1 <- g(x) / f1(x)
gf2 <- g(x) / f2(x)
plot(x, gf1, type = 'l', xlim = c(1,5), lty = 2, col = 2, lwd = 2)
lines(x, gf2, lty = 3, col = 3, lwd = 2)
legend("topright", c(expression(g(x)/f1(x)), expression(g(x)/f2(x))),
       inset = 0.02, lty = 2:3, col = 2:3, lwd = 2)

## ----comment=''---------------------------------------------------------------
set.seed(22034)
m <- 10000
x <- rexp(m, 1) + 1
gf1 <- g(x) / f1(x)
x <- abs(rnorm(m)) + 1
gf2 <- g(x) / f2(x)
theta.hat = c(mean(gf1), mean(gf2))
var.hat = c(var(gf1), var(gf2)) / m
rbind(theta.hat, var.hat)

## ----comment=''---------------------------------------------------------------
set.seed(22034)
m <- 10000
g <- function(x) {exp(-x) / (1 + x^2) * (x>0) * (x<1)}
f <- function(x) {exp(-x) / (1 - exp(-1)) * (x>0) * (x<1)}
u <- runif(m)
x <- - log(1 - u * (1 - exp(-1))) # inverse transformation method
gf <- g(x) / f(x)
var1 <- var(gf)
c(mean(gf), var1)

## ----comment=''---------------------------------------------------------------
set.seed(22034)
M <- 10000
k <- 5
m <- M/k
theta.hat <- var.hat <- numeric(k)
g <- function(x) {exp(-x) / (1 + x^2) * (x>0) * (x<1)}
f <- function(x) {k * exp(-x) / (1 - exp(-1)) * (x>0) * (x<1)}
for (j in 1:k) {
  u <- runif(m, (j-1)/k, j/k)
  x <- -log(1 - u * (1 - exp(-1)))
  gf <- g(x) / f(x)
  theta.hat[j] <- mean(gf)
  var.hat[j] <- var(gf)
}
var2 <- sum(var.hat) # the variance of theta is sum of the variance of each stratum
c(sum(theta.hat), var2)

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----comment=''---------------------------------------------------------------
# clear memory, set seed, and initialize sample size
rm(list = ls())
set.seed(22034)
n <- 100

# construct confidence interval
CI_t <- function(n) {
  # generate data
  x <- rlnorm(n)
  y <- log(x)
  
  # calculate the estimation of parameter
  mu.hat <- mean(y)
  sigma.hat <- sd(y)
  
  # calculate confidence interval (note we construct CI with t distribution)
  result <- mu.hat + sigma.hat / sqrt(n) * qt(c(0.025, 0.975), n - 1) 
  return(result)
}

# replicate the function
m <- 10000
MC.CI <- matrix(nrow = m, ncol = 2)
for (i in 1:m) { MC.CI[i,] <- CI_t(n) }
mean(MC.CI[,1] < 0 & MC.CI[,2] > 0)

## -----------------------------------------------------------------------------
# clear memory, and set seed
rm(list = ls())
set.seed(22034)

# Count Five Test
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outX <- sum(X > max(Y)) + sum(X < min(Y)) # extreme obsevations of X to Y
  outY <- sum(Y > max(X)) + sum(Y < min(X)) # extreme obsevations of Y to X
  return(as.integer(max(c(outX, outY)) > 5)) # 1 for rejecting the null hypothesis
}

## -----------------------------------------------------------------------------
# intitialize
n_test <- c(10, 20, 50, 100, 200, 500) # number of sample size
m <- 1000 # numebr of iterations
sigmaX <- 1
sigmaY <- 1.5

# conduct the experiment
result <- matrix(nrow = 6, ncol = 3)
for (n in n_test) {
  result_n <- matrix(nrow = m, ncol = 2)
  for (i in 1:m) {
    # generate data 
    x <- rnorm(n, 0, sigmaX)
    y <- rnorm(n, 0, sigmaY)
    
    # conduct the two test
    count5test.res <- count5test(x, y)
    Ftest.res <- as.integer(var.test(x, y)$p.value <= 0.055)
    
    # save the result
    result_n[i,] <- c(count5test.res, Ftest.res)
  }
  result[which(n_test == n),] <- c(n, apply(result_n, 2, mean))
}

## ----comment=''---------------------------------------------------------------
print(result)

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----comment=''---------------------------------------------------------------
# clear memory, import package, and attach dataset
rm(list = ls())
library(boot)
xsample <- aircondit$hours

# Bootstrap
set.seed(22034)
B <- 10000
lambda <- 1 / mean(xsample)
lambda.star <- numeric(B)
for (b in 1:B) {
  xstar <- sample(xsample, replace = TRUE)
  lambda.star[b] <- 1 / mean(xstar)
}

# reuslt
round(c(lambda = lambda,
        lambda.star = mean(lambda.star),
        bias = mean(lambda.star) - lambda,
        se.boot = sd(lambda.star)), 4)

## ----comment=''---------------------------------------------------------------
# Bootstrap with the package boot
set.seed(22034)
hazard_rate <- function(x, i) { 1 / mean(x[i]) }
boot.res <- boot(xsample, statistic = hazard_rate, R = 10000)

# result_standard
round(c(original = boot.res$t0,
        bias = mean(boot.res$t) - boot.res$t0,
        se.boot = sd(boot.res$t)), 4)

## ----comment=''---------------------------------------------------------------
# Bootstrap
set.seed(22034)
B <- 10000
meantime <- mean(xsample)
meantime.star <- numeric(B)
for (b in 1:B) {
  xstar <- sample(xsample, replace = TRUE)
  meantime.star[b] <- mean(xstar)
}

# reuslt
round(c(meantime = meantime,
        meantime.star = mean(meantime.star),
        bias = mean(meantime.star) - meantime,
        se.boot = sd(meantime.star)), 4)

## ----comment=''---------------------------------------------------------------
CI.std <- meantime + qnorm(c(0.025, 0.975)) * sd(meantime.star)

## ----comment=''---------------------------------------------------------------
CI.basic <- 2 * meantime - unname(quantile(meantime.star, c(0.975, 0.025)))

## ----comment=''---------------------------------------------------------------
CI.perc <- unname(quantile(meantime.star, c(0.025, 0.975)))

## ----comment=''---------------------------------------------------------------
# the bias correction factor
z0 <- qnorm(sum(meantime.star < meantime) / B)

# the acceleratiom factor (jackknife est.)
meantime.jack <- numeric(B)
for (i in 1:B) {
  meantime.jack[i] <- mean(xsample[-i])
}
L <- mean(meantime.jack) - meantime.jack
a <- sum(L^3) / (6 * sum(L^2)^1.5)

# BCa conf. limits
alpha <- c(0.025, 0.975)
zalpha <- qnorm(alpha)
adj.alpha <- pnorm(z0 + (z0 + zalpha) / (1 - a * (z0 + zalpha)))
CI.BCa <- unname(quantile(meantime.star, adj.alpha))

## ----comment=''---------------------------------------------------------------
rbind(CI.std, CI.basic, CI.perc, CI.BCa)

## -----------------------------------------------------------------------------
# Bootstrap with the package boot
set.seed(22034)
mean_time <- function(x, i) { mean(x[i]) }
boot.res <- boot(xsample, statistic = mean_time, R = 10000)
boot.ci(boot.res, type = c("norm", "basic", "perc", "bca"))

## ----fig.height=4, fig.width=6------------------------------------------------
hist(boot.res$t, breaks = "scott", prob = TRUE)
abline(v = boot.res$t0, col = 'red', lwd = 2)

## ----comment=''---------------------------------------------------------------
# clear memory, import package, and set seed
rm(list = ls())
library(boot)
set.seed(22034)

# Bootstrap
m <- 1000
mu <- 1
boot.mean <- function(x, i) { mean(x[i]) }
CI.norm <- CI.basic <- CI.perc <- matrix(nrow = m, ncol = 2)
for (i in 1:m) {
  xsample <- rnorm(100, mu, 1)
  res <- boot(data = xsample, statistic = boot.mean, R = 1000)
  ci <- boot.ci(res, type = c("norm", "basic", "perc"))
  CI.norm[i,] <- ci$norm[2:3]
  CI.basic[i,] <- ci$basic[4:5]
  CI.perc[i,] <- ci$perc[4:5]
}

# result
coverage_rate <- c(norm = mean(CI.norm[, 1] <= mu & CI.norm[, 2] >= mu),
                   basic = mean(CI.basic[, 1] <= mu & CI.basic[, 2] >= mu),
                   perc = mean(CI.perc[, 1] <= mu & CI.perc[, 2] >= mu))
miss_left <- c(mean(CI.norm[, 1] > mu), mean(CI.basic[, 1] > mu), mean(CI.perc[, 1] > mu))
miss_right <- c(mean(CI.norm[, 2] < mu), mean(CI.basic[, 2] < mu), mean(CI.perc[, 2] < mu))
rbind(coverage_rate, miss_left, miss_right)

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----comment=''---------------------------------------------------------------
# clear memory, import package, and set seed
rm(list = ls())
library(bootstrap)
set.seed(22034)

# attach the dataset
attach(scor)
x <- data.frame(scor)
detach(scor)
n <- nrow(x)

# jackknife
theta.jack <- numeric(n)
lambda <- eigen(cov(x))$values
theta.hat <- max(lambda) / sum(lambda)
for (i in 1:n) {
  y <- x[-i,]
  lambda.jack <- eigen(cov(y))$values
  theta.jack[i] <- max(lambda.jack) / sum(lambda.jack)
}

# result
theta.jack.mean <- mean(theta.jack)
bias <- (n - 1) * (theta.jack.mean - theta.hat)
se.jack <- sqrt((n - 1) * mean((theta.jack - theta.jack.mean)^2))
c(theta.hat = theta.hat, theta.jack = theta.jack.mean, bias = bias, se.jack = se.jack)

## ----comment='', message=F----------------------------------------------------
# clear memory, set seed, and import package
rm(list = ls())
set.seed(22034)
library(DAAG)

# attach the data set
attach(ironslag)
n <- length(magnetic)

# leave-two-out
e1 <- e2 <- e3 <- e4 <- numeric(choose(n, 2))
ij <- 1
for (i in 1:(n-1)) 
  for (j in (i+1):n) {
    # combine the two indices together
    k <- c(i, j)
    y <- magnetic[-k]
    x <- chemical[-k]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
    e1[ij] <- mean((magnetic[k] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
    e2[ij] <- mean((magnetic[k] - yhat2)^2)
    
    J3 <- lm(log(y) ~ x)
    yhat3 <- exp(J3$coef[1] + J3$coef[2] * chemical[k])
    e3[ij] <- mean((magnetic[k] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    yhat4 <- exp(J4$coef[1] + J4$coef[2] * log(chemical[k]))
    e4[ij] <- mean((magnetic[k] - yhat4)^2)
    
    ij <- ij + 1
  }

# result
c(Linear = mean(e1), Quad = mean(e2), Exp = mean(e3), LogLog = mean(e4))

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
L2 <- lm(magnetic ~ chemical + I(chemical^2))
detach(ironslag)
plot(L2)

## ----message=F----------------------------------------------------------------
# clear memory, set seed, and import package
rm(list = ls())
set.seed(22034)
library(MASS)

# spearman rank correlation test
spearman.test <- function(x, y, R = 49999) {
  S0 <- cor.test(x, y, method = "spearman")$estimate
  S <- numeric(R)
  for (i in 1:R) {
    k <- sample(1:length(x))
    S[i] <- cor.test(x, y[k], method = "spearman")$estimate
  }
  p.value <- mean(c(S0, S) >= S0)
  return(c(statistic = S0, p.value = p.value))
}

## -----------------------------------------------------------------------------
mu <- c(0, 1)
sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
n <- 50
R <- 49999
x <- mvrnorm(n, mu, sigma)
spearman.test(x[, 1], x[, 2], R = R)
cor.test(x[, 1], x[, 2], method = "spearman")

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22034)

rl.metropolis <- function(sigma, x0, N) {
  # sigma: sd of proposal distribution N(xt,sigma^2)
  # x0: initial value
  # N: length of chain
  
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0  # to calculate acceptance rate
  for (t in 2:N) {
    y <- rnorm(1, x[t-1], sigma)
    if (u[t] <= exp(abs(x[t-1]) - abs(y))) { x[t] <- y; k <- k + 1 }
    else { x[t] <- x[t-1] }
  }
  return(list(mc = x, acc.prob = k / N))
}

N <- 10000
b <- 1000
k <- 4
sigma <- c(0.5, 1, 4, 16)
x0 <- c(-5, -2, 2, 5)
X <- matrix(nrow = k, ncol = N)
acc.prob <- numeric(k)
for (i in 1:k) {
  rl <- rl.metropolis(sigma[i], x0[i], N)
  X[i, ] <- rl$mc
  acc.prob[i] <- rl$acc.prob
}
acc.prob

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
for (i in 1:k) {
  plot(X[i,], type = "l", xlab = bquote(sigma == .(sigma[i])),
       ylab = "X", ylim = range(X[i,]))
}

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2, 2))
x <- seq(-6, 6, 0.01)
fx <- exp(-abs(x)) / 2
for (i in 1:k) {
  hist(X[i, -(1:b)], breaks = "Scott", freq = FALSE, main = "",
       xlab = bquote(sigma == .(sigma[i])), xlim = c(-6, 6), ylim = c(0, 0.5),)
  lines(x, fx, col = 2, lty = 2)
}

## -----------------------------------------------------------------------------
z <- rexp(100, 1)
z <- c(-rev(z), z) # generate laplace random numbers
p <- c(0.05, seq(0.1, 0.9, 0.1), 0.95)
Q <- quantile(z, p)
mc <- X[, -(1:b)]
Qmc <- apply(mc, 1, function(x) quantile(x, p))
QQ <- data.frame(round(cbind(Q, Qmc), 3))
names(QQ) <- c('True', 'sigma=0.5', 'sigma=1', 'sigma=4', 'sigma=16')
knitr::kable(QQ)

## ----fig.height=4, fig.width=6------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# ergodic mean plot
phi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(phi)) {
  phi[i,] <- phi[i,] / (1:ncol(phi))
}
for (i in 1:k) {
  if (i == 1) {
    plot((b+1):N, phi[i, (b+1):N], ylim = c(-0.5, 0.5),
         type = "l", xlab = 'Index', ylab = bquote(phi))
  } else { lines(phi[i, (b+1):N], col = i) }
}

# plot of R_hat
rhat <- rep(0, N)
for (j in (b+1):N) {
  rhat[j] <- Gelman.Rubin(phi[, 1:j])
}
plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R")
abline(h = 1.2, lty = 2)

## ----fig.height=4, fig.width=6------------------------------------------------
# clear memory and set seed
#rm(list = ls())
set.seed(22034)

rbn.metropolis <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XY <- rbn.metropolis(mu, sigma, rho, mu, N)
X <- XY$X[-(1:b)]; Y <- XY$Y[-(1:b)]
plot(X, Y, xlab = bquote(X[t]), ylab = bquote(Y[t]),
     main = "", cex = 0.5, ylim = range(Y))
cov(cbind(X, Y))

## ----fig.height=4, fig.width=9------------------------------------------------
k <- 4
x0 <- matrix(c(2,2,-2,-2,4,-4,-4,4), nrow = 2, ncol = k)
Xmc <- Ymc <- XYmc <- matrix(0, nrow = k, ncol = N)
for (i in 1:k) {
  XY <- rbn.metropolis(mu, sigma, rho, x0[,i], N)
  Xmc[i,] <- XY$X; Ymc[i,] <- XY$Y
  XYmc[i,] <- Xmc[i,] * Ymc[i,]
}

# ergodic mean plot
cal_phi <- function(X) {
  phi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(phi)) {
    phi[i,] <- phi[i,] / (1:ncol(phi))
  }
  return(phi)
}
phiX <- cal_phi(Xmc)
phiY <- cal_phi(Ymc)
phiXY <- cal_phi(XYmc)

plot_erg_mean <- function(phi, rg) {
  for (i in 1:k) {
    if (i == 1) {
      plot((b+1):N, phi[i, (b+1):N], type = "l", ylim = rg,
           xlab = "Index", ylab = bquote(phi))
    }
    else { lines(phi[i, (b+1):N], col = i) }
  }
}
par(mfrow = c(1, 3))
plot_erg_mean(phiX, rg = c(-0.5, 0.5))
plot_erg_mean(phiY, rg = c(-0.5, 0.5))
plot_erg_mean(phiXY, rg = c(0.7, 1.1))

## ----fig.height=4, fig.width=9------------------------------------------------
Gelman.Rubin <- function(phi) {
  phi <- as.matrix(phi)
  k <- nrow(phi); n <- ncol(phi)
  phi.means <- rowMeans(phi)
  B <- n * var(phi.means)
  phi.w <- apply(phi, 1, var)
  W <- mean(phi.w)
  v.hat <- W * (n - 1) / n + B / n
  r.hat <- v.hat / W
  return(r.hat)
}

# plot of R_hat
plot_R_hat <- function(phi) {
  rhat <- rep(0, N)
  for (j in (b+1):N) {
    rhat[j] <- Gelman.Rubin(phi[, 1:j])
  }
  plot(rhat[(b+1):N], type = "l", xlab = "", ylab = "R", ylim = c(1, 1.25))
  abline(h = 1.2, lty = 2)
}
par(mfrow = c(1, 3))
plot_R_hat(phiX)
plot_R_hat(phiY)
plot_R_hat(phiXY)

## ----comment = ''-------------------------------------------------------------
lm.fit <- lm(Y ~ X)
summary(lm.fit)

## ----fig.height=5, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
e <- lm.fit$residuals
qx <- seq(-2, 2, 0.01)
hist(e, breaks = "Scott", freq = FALSE, main = "", xlim = c(-2, 2), ylim = c(0, 1))
lines(qx, dnorm(qx, 0, sqrt(0.19)), col = 2, lwd = 1.5)
qqnorm(e)
qqline(e, col = 2, lwd = 2, lty = 2)

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## ----message=FALSE, warning=FALSE---------------------------------------------
# clear memory, import package and set seed
rm(list = ls())
library(mediation)
set.seed(22034)

# generate data
gen_data <- function(n, aM, aY, alpha, beta, gamma) {
  X <- rnorm(n, 0, 1)
  M <- aM + alpha * X + rnorm(n)
  Y <- aY + beta * M + gamma * X + rnorm(n)
  return(data.frame(X = X, M = M, Y = Y))
}
n <- 100
aM <- 1; aY <- 2
alpha <- beta <- gamma <- 1
dataset <- gen_data(n, aM, aY, alpha[1], beta[1], gamma)

## -----------------------------------------------------------------------------
# mediation effect estimate
model.Y <- lm(Y ~ X + M, data = dataset)
model.M <- lm(M ~ X, data = dataset)
res <- mediate(model.M, model.Y, treat = "X", mediator = "M", sims = 10)
summary(res)

## -----------------------------------------------------------------------------
alphabeta.hat <- res$d0
se.hat <- sd(res$d0.sims)
T0 <- alphabeta.hat / se.hat
T0

## -----------------------------------------------------------------------------
med.perm.test <- function(alpha, beta, R = 499) {
  dataset <- gen_data(n, aM, aY, alpha, beta, gamma)
  
  # origin full models
  model.Y <- lm(Y ~ X + M, data = dataset)
  model.M <- lm(M ~ X, data = dataset)
  res <- mediate(model.M, model.Y, treat = "X", mediator = "M", sims = 10)
  p.value.mediate <- res$d0.p
  alphabeta.hat <- res$d0
  T0 <- res$d0 / sd(res$d0.sims)
  
  # reduced models
  model.Yr <- lm(Y ~ X, data = dataset)
  Y.hat <- model.Yr$fitted.values
  e.Yr <- model.Yr$residuals
  model.Mr <- lm(M ~ 1, data = dataset)
  M.hat <- model.Mr$fitted.values
  e.Mr <- model.Mr$residuals
  
  # permuation
  Tstar <- numeric(R)
  for (i in 1:R) {
    # calculate new data
    e.Ystar <- e.Yr[sample(1:n, replace = FALSE)]
    e.Mstar <- e.Mr[sample(1:n, replace = FALSE)]
    Ystar <- Y.hat + e.Ystar
    Mstar <- M.hat + e.Mstar
    
    # new full models
    datastar <- data.frame(X = dataset$X, Ystar = Ystar, Mstar = Mstar)
    model.Ystar <- lm(Ystar ~ X + Mstar, data = datastar)
    model.Mstar <- lm(Mstar ~ X, data = datastar)
    res <- mediate(model.Mstar, model.Ystar, treat = "X", mediator = "Mstar", sims = 10)
    Tstar[i] <- res$d0 / sd(res$d0.sims)
  }
  p.value.perm <- mean(c(abs(Tstar), abs(T0)) >= abs(T0))
  return(c(alphabeta.hat = alphabeta.hat, statistic = T0,
           p.value.perm = p.value.perm, p.value.mediate = p.value.mediate))
}

## -----------------------------------------------------------------------------
med.perm.test(alpha = 1, beta = 1, R = 999)
med.perm.test(alpha = 1, beta = 0, R = 999)
med.perm.test(alpha = 0, beta = 1, R = 999)
med.perm.test(alpha = 0, beta = 0, R = 999)

## -----------------------------------------------------------------------------
# clear memory, and set seed
rm(list = ls())
set.seed(22034)

# make a function
get_alpha <- function(N, b, f0) {
  x <- cbind(X1 = rpois(N, 1), X2 = rexp(N, 1),
             X3 = sample(0:1, N, replace = TRUE))
  g <- function(alpha) {
    p <- 1 / (1 + exp(-alpha - x %*% b))
    return(mean(p) - f0)
  }
  res <- uniroot(g, c(-20, 0))
  return(res$root)
}

## ----fig.height=4, fig.width=6------------------------------------------------
# calculate and plot
f0 <- c(0.1, 0.01, 0.001, 0.0001)
alpha <- numeric(4)
for (i in 1:4) {
  alpha[i] <- get_alpha(1e6, c(0, 1, -1), f0[i])
}
plot(f0, alpha)

## ----fig.height=4, fig.width=6------------------------------------------------
# calculate and plot again
f0 <- seq(0.0001, 0.1, 0.001)
alpha <- numeric(length(f0))
for (i in 1:length(f0)) {
  alpha[i] <- get_alpha(1e4, c(2, 1, -1), f0[i])
}
plot(f0, alpha, type = "l")

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

# create function to calc s(lambda)
s.lambda <- function(lambda) {
  u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
  v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
  return(sum((u * exp(-lambda * u) - v * exp(-lambda * v)) /
               (exp(-lambda * u) - exp(-lambda * v))))
}

# using `uniroor` function
lambda.mle <- uniroot(s.lambda, c(0, 10))$root
round(lambda.mle, 5)

## -----------------------------------------------------------------------------
# create function to calc s'(lambda)
ds.lambda <- function(lambda) {
  u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
  v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
  return(sum((u - v)^2 * exp(-lambda * (u + v)) / 
               (exp(-lambda * u) - exp(-lambda * v))^2))
}

# using newton's method
lambda0 <- 0.1
while (TRUE) {
  lambda1 <- lambda0 - s.lambda(lambda0) / ds.lambda(lambda0)
  if (abs(lambda1 - lambda0) > 1e-6) { lambda0 <- lambda1 } else { break }
}
round(lambda1, 5)

## -----------------------------------------------------------------------------
x <- list(c(1, 2), list(3, 4))
str(x)
str(unlist(x))
str(as.vector(x))

## -----------------------------------------------------------------------------
dim(c(1, 2, 3)) # atomic vector
dim(list(1, 2, list(3))) # list

## -----------------------------------------------------------------------------
x <- matrix(1:6, nrow = 2, ncol = 3)
c(is.matrix(x), is.array(x))

y <- array(1:12, c(2, 3, 2))
c(is.matrix(y), is.array(y))

## -----------------------------------------------------------------------------
x <- data.frame(V1 = c(1, 2, 3),
                V2 = c("a", "b", "c"),
                V3 = c(TRUE, FALSE, FALSE),
                row.names = c("X1", "X2", "X3"))
x
attributes(x)
dim(x)

## -----------------------------------------------------------------------------
x <- data.frame(
  V1 = c(1L, 2L),
  V2 = c(FALSE, TRUE),
  V3 = c("a", "b")
)
as.matrix(x)

y <- data.frame(
  V1 = c(1.5, 2.0),
  V2 = c(FALSE, TRUE)
)
as.matrix(y)

## -----------------------------------------------------------------------------
# 0 rows
x <- data.frame(V1 = numeric())
c(nrow(x), ncol(x))

# 0 columns
y <- data.frame(row.names = "X1")
c(nrow(y), ncol(y))

# 0 rows, 0 columns
z <- data.frame()
c(nrow(z), ncol(z))

## -----------------------------------------------------------------------------
# clear memory
rm(list = ls())

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
x <- data.frame(x1 = c(1.5, 2.5, 3.5, 4.5), x2 = rnorm(4, 4, 4))
str(x)

## -----------------------------------------------------------------------------
res1 <- data.frame(lapply(x, scale01))
res1

## -----------------------------------------------------------------------------
# add a non-numeric column
x$x3 = c(rep("A", 2), rep("B", 2))

res2 <- data.frame(lapply(x, function(x) if (is.numeric(x)) scale01(x) else x))
res2

## -----------------------------------------------------------------------------
rm(list = ls())
x <- data.frame(x1 = c(0.6, 1.3, 7.6, 2.4), x2 = rnorm(4, 2, 2))
str(x)

## -----------------------------------------------------------------------------
res1 <- vapply(x, sd, 1)
res1

## -----------------------------------------------------------------------------
# add a non-numeric column
x$x3 = c(rep("A", 2), rep("B", 2))

res2 <- vapply(x[vapply(x, is.numeric, TRUE)], sd, 1)
res2

## -----------------------------------------------------------------------------
# clear memory and set seed
rm(list = ls())
set.seed(22034)

gibbsR <- function(mu, sigma, rho, initial, N) {
  # mu, sigma, rho: parameter of bivariate normal distribution.
  # initial: initial value
  # N: length of chain
  
  X <- Y <- numeric(N)
  s <- sqrt(1 - rho^2) * sigma
  X[1] <- initial[1]; Y[1] <- initial[2]
  for (i in 2:N) {
    y <- Y[i-1]
    m1 <- mu[1] + rho * (y - mu[2]) * sigma[1] / sigma[2]
    X[i] <- rnorm(1, m1, s[1])
    x <- X[i]
    m2 <- mu[2] + rho * (x - mu[1]) * sigma[2] / sigma[1]
    Y[i] <- rnorm(1, m2, s[2])
  }
  return(list(X = X, Y = Y))
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction("NumericMatrix gibbsC(NumericVector mu, NumericVector sigma, double rho,
                     NumericVector initial, int N) {
  // mu, sigma, rho: parameter of bivariate normal distribution.
  // initial: initial value
  // N: length of chain
  
  NumericMatrix XY(N, 2);
  double x, y, m1, m2;
  XY(0, 0) = initial[0];
  XY(0, 1) = initial[1];
  for(int i = 1; i < N; i++) {
    y = XY(i - 1, 1);
    m1 = mu[0] + rho * (y - mu[1]) * sigma[0] / sigma[1];
    XY(i, 0) = rnorm(1, m1, sqrt(1 - rho * rho))[0] * sigma[0];
    x = XY(i, 0);
    m2 = mu[1] + rho * (x - mu[0]) * sigma[1] / sigma[0];
    XY(i, 1) = rnorm(1, m2, sqrt(1 - rho * rho))[0] * sigma[1];
  }
  return(XY);
}")

## ----fig.height=4, fig.width=8------------------------------------------------
# generate chains
N <- 10000
b <- 1000
rho <- 0.9
mu <- c(0, 0)
sigma <- c(1, 1)
XYR <- gibbsR(mu, sigma, rho, mu, N)
XR <- XYR$X[-(1:b)]; YR <- XYR$Y[-(1:b)]
#sourceCpp('gibbsC.cpp')
XYC <- gibbsC(mu, sigma, rho, mu, N)
XC <- XYC[-(1:b), 1]; YC <- XYC[-(1:b), 2]

par(mfrow = c(1, 2))
qqplot(XR, XC, plot.it = TRUE)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
qqplot(YR, YC, plot.it = TRUE)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)

## ----fig.height=4, fig.width=8------------------------------------------------
par(mfrow = c(1, 2))
plot(XR, YR, cex = 0.5)
plot(XC, YC, cex = 0.5)

## -----------------------------------------------------------------------------
# import package
library(microbenchmark)
ts <- microbenchmark(gibbsR = gibbsR(mu, sigma, rho, mu, N),
                     gibbsC = gibbsC(mu, sigma, rho, mu, N))
summary(ts)[, c(1, 3, 5, 6)]

