# Consider a terminal node which has 12 observations in it
# There are 2 continents (using e.g. the gap minder data)
# with 2 countries in the first continent and 3 countries in the second
# continent

# The data looks like this:
df <- data.frame(
  obs = 1:12,
  R = c(0.942, -0.751, -0.993, -2.15, -0.906, -1.94, -0.24, -0.741,
        -0.967, 1.996, 0.134, 0.325),
  continent = factor(c(rep(1, 6), rep(2, 6))),
  country = factor(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5)), 
  x1 = runif(12)
)

# Need a matrix for allocation of nesting
# So first 2 countries are in continent 1
# Next 3 countries are in continent 2
nesting <- c(1, 1, 2, 2, 2)

# We can now form the matrices M1 and M2
M1 <- model.matrix(~ df$country - 1)
M2 <- model.matrix(~ df$continent - 1)

# Now suppose we have some parameter values
tau <- 1
tau_mu <- 1
tau_phi <- 1
tau_gamma <- 1
P <- 1

# We need these matrices:
n <- nrow(df)
n_phi <- nlevels(df$continent)
n_gamma <- nlevels(df$country)
M1M1T <- tcrossprod(M1)
M2M2T <- tcrossprod(M2)
Delta_R <- diag(n)/tau + M1M1T/(P * tau_gamma)
Phi_R <- Delta_R + M2M2T/(P * tau_phi)

# Update mu
prec <- sum(solve(Phi_R)) + tau_mu
mean <- sum(solve(Phi_R, df$R)) / prec
mu <- rnorm(1, mean = mean, sd = sqrt( 1 / prec ))

# Update phi
prec <- t(M2)%*%solve(Delta_R, M2) + (P*tau_phi) * diag(n_phi)
mean <- solve(prec, (t(M2)%*%solve(Delta_R, df$R) + (P*tau_phi * mu)))
phi <- mvnfast::rmvn(1, mu = mean,
                     sigma = solve(prec))

# Update gamma
prec <- tau * crossprod(M1) + tau_gamma * P * diag(n_gamma)
mean <- solve(prec, tau * t(M1) %*% df$R + tau_gamma * P * phi[nesting])
gamma <- mvnfast::rmvn(1, mu = mean,
                       sigma = solve(prec))


