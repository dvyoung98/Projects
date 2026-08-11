#install needed package
install.packages("stochvol")
library(stochvol)
citation(stochvol)

#set seed and visualize data
set.seed(123)
data(exrates)
ret <- logret(exrates$USD, demean = TRUE)
par(mfrow = c(2, 1), mar = c(1.9, 1.9, 1.9, 0.5), mgp = c(2, 0.6, 0))
plot(exrates$date, exrates$USD, type = "l",
     main = "Price of 1 EUR in USD")
plot(exrates$date[-1], ret, type = "l", main = "Demeaned log returns")

#simulate and plot time series
sim <- svsim(500, mu = -9, phi = 0.99, sigma = 0.1)
par(mfrow = c(2, 1))
plot(sim)


#playing around, running sampler
ret <- ret[1:200]
res <- svsample(ret, priormu = c(-10, 1), priorphi = c(20, 1.1),
                priorsigma = 0.1, thin = 10)

summary(res, showlatent = FALSE)


res <- updatesummary(res, quantiles = c(0.01, 0.1, 0.5, 0.9, 0.99))
volplot(res, forecast = 100, dates = exrates$date[seq_along(ret)])


par(mfrow = c(1, 3))
paradensplot(res, showobs = FALSE)



install.packages("sandwich")
library(sandwich)




#Bayesian normal linear model with homoskedastic errors
#Set seed
set.seed(123456)
n <- 200
beta.true <- c(0.1, 0.5)
sigma.true <- 0.01
X <- matrix(c(rep(1, n), rnorm(n, sd = sigma.true)), nrow = n)
y <- rnorm(n, X %*% beta.true, sigma.true)

#set parameters
burnin <- 100
draws <- 1000
b0 <- matrix(c(0, 0), nrow = ncol(X))
B0inv <- diag(c(10^-10, 10^-10))
c0 <- 0.001
C0 <- 0.001

#calculate outside MCMC main loop
p <- ncol(X)
preCov <- solve(crossprod(X) + B0inv)
preMean <- preCov %*% (crossprod(X, y) + B0inv %*% b0)
preDf <- c0 + n/2 + p/2

#assign storage space
draws1 <- matrix(NA_real_, nrow = draws, ncol = p + 1)
colnames(draws1) <- c(paste("beta", 0:(p-1), sep = "_"), "sigma")
sigma2draw <- 1

#Run main sampler
for (i in -(burnin-1):draws) {
  betadraw <- as.numeric(mvtnorm::rmvnorm(1, preMean,
                                          sigma2draw * preCov))
  tmp <- C0 + 0.5 * (crossprod(y - X %*% betadraw) + 
                       crossprod((betadraw - b0), B0inv) %*% (betadraw - b0))
  sigma2draw <- 1 / rgamma(1, preDf, rate = tmp)
  if (i > 0) draws1[i, ] <- c(betadraw, sqrt(sigma2draw))
}

#obtain point estimates for parameters
colMeans(draws1)

#plot
plot(coda::mcmc(draws1))




#Bayesian normal linear model with SV errors
#simulate
mu.true <- log(sigma.true^2)
phi.true <- 0.97
vv.true <- 0.3
simresid <- svsim(n, mu = mu.true, phi = phi.true, sigma = vv.true)
y <- X %*% beta.true + simresid$y

#parameters
draws <- 5000
burnin <- 500
thinning <- 10
priors <- specify_priors(
  mu = sv_normal(-10, 2),
  phi = sv_beta(20, 1.5),
  sigma2 = sv_gamma(0.5, 0.5)
  )

#assign storage space
draws2 <- matrix(NA_real_, nrow = floor(draws / thinning),
                 ncol = 3 + n + p)
colnames(draws2) <- c("mu", "phi", "sigma",
                      paste("beta", 0:(p-1), sep = "_"), paste("h", 1:n, sep = "_"))
betadraw <- c(0, 0)
paradraw <- list(mu = -10, phi = 0.9, sigma = 0.2)
latentdraw <- rep(-10, n)
paranames <- names(paradraw)

#Run main sampler
for (i in -(burnin-1):draws) {
  ytilde <- y - X %*% betadraw
  svdraw <- svsample_fast_cpp(ytilde, startpara = paradraw,
                              startlatent = latentdraw, priorspec = priors)
  paradraw <- svdraw$para
  latentdraw <- drop(svdraw$latent)
  normalizer <- as.numeric(exp(-latentdraw / 2))
  Xnew <- X * normalizer
  ynew <- y * normalizer
  Sigma <- solve(crossprod(Xnew) + B0inv)
  mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
  betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))
  if (i > 0 && i %% thinning == 0) {
    draws2[i/thinning, 1:3] <- drop(paradraw)[paranames]
    draws2[i/thinning, 4:5] <- betadraw
    draws2[i/thinning, 6:(n+5)] <- latentdraw
  }
}

#plot
plot(coda::mcmc(draws2[, 4:8]))
