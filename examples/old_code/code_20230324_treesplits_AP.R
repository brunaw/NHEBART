# Simulate data for a nested model then fit it using the simple maths
# in code_20230323.R
# This version has a simple tree in it too

# Load in used packages
library(tidyverse)
library(mvnfast)
library(progress)
devtools::load_all(".")
# Simulate the data -------------------------------------------------------

set.seed(100)
n <- 500
# So the below is like the number of continents
n_groups <- 4
# And the below is the number of countries in each continent
n_groups_nested <- c(3, 2, 5, 3)
group <- sample(1:n_groups, size = n, replace = TRUE)
sub_group <- rep(NA, length = n)
for(i in 1:length(group)) {
  sub_group[i] <- sample(1:n_groups_nested[group[i]],
                         size = 1)
}

# Might want to fiddle around with these values a bit
tau_mu <- 3
tau_phi <- 3
tau_gamma <- 2
tau <- 1

# Now set up the tree
n_terminal_nodes <- 2
x <- runif(n)
# The split will be at x < 0.5

# Simulate terminal node parameter
true_mu <- rnorm(n_terminal_nodes, 0, 1/sqrt(tau_mu))

# Simulate grouping parameter
true_phi <- vector('list', length = n_terminal_nodes)
for (i in 1:n_terminal_nodes) {
  true_phi[[i]] <- rnorm(n_groups, true_mu, 1/sqrt(tau_phi))
}

# Simulate sub-grouping parameters
true_gamma <- vector('list', length = n_terminal_nodes)
for (i in 1:n_terminal_nodes) {
  true_gamma[[i]] <- vector('list', length = n_groups)
  for(j in 1:length(true_gamma[[i]])) {
    true_gamma[[i]][[j]] <- rnorm(n_groups_nested[j], true_phi[[i]][j], 1/sqrt(tau_gamma))
  }
}

# Introduce the tree
# For phi the order of the indices is terminal node, group
# So phi[[1]][3] is the 1st terminal node, 3rd group
# For gamma it's the same but with sub-group added on
# lambda[[1]][[2]][3] is the 1st terminal node, 2nd group,
# 3rd sub-group

tree_fun <- function(x, group, subgroup) {
  ans <- if(x < 0.5 ) {
    true_gamma[[1]][[group]][subgroup]
  } else {
    true_gamma[[2]][[group]][subgroup]
  }
  return(ans)
}

# Now simulate the data
y <- rep(NA, length = n)
for (i in 1:n) {
  y[i] <- rnorm(1, mean =
                  tree_fun(x[i], group[i], sub_group[i]),
                sd = 1 / sqrt(tau))
}

# Get the data into the nice format from yesterday ------------------------

# The data looks like this:
df <- data.frame(
  obs = 1:n,
  R = y,
  x = x,
  continent = factor(group),
  country = factor(paste(group,sub_group,sep = "_"))
)

# Need a matrix for allocation of nesting
nesting <- rep(1:n_groups, n_groups_nested)

# Now use code from 23/3/23 to fit the model ------------------------------

# Now run the MCMC
n_iter <- 500

# Store everything for each terminal node
mu_store <- vector('list', length = n_terminal_nodes)
phi_store <- vector('list', length = n_terminal_nodes)
gamma_store <- vector('list', length = n_terminal_nodes)
for (i in 1:length(mu_store)) {
  mu_store[[i]] <- rep(NA, length = n_iter)
  phi_store[[i]] <- matrix(NA, nrow = n_iter, ncol = n_groups)
  gamma_store[[i]] <- matrix(NA, nrow = n_iter, ncol = length(nesting))
}

pb = progress::progress_bar$new(
  format = "Progress [:bar]  :current/:total (:percent)",
  width = 60, total = n_iter)

# Pick out the two data sets
which_left <- which(x < 0.5)
which_right <- which(x >= 0.5)
M1 <- model.matrix(~ df$country - 1) # Note the column names are a bit confusing
M2 <- model.matrix(~ df$continent - 1)
P <- 1 # Number of trees
formula = R ~ x
group_variable = "continent"
subgroup_variable = "country"
num_trees <- 1

tau         = 1
tau_phi     = 3 
tau_mu      = 3
tau_gamma   = 2
data = df


for (i in 1:n_iter) {
  pb$tick()

  # Do it for each terminal node
  for (j in 1:n_terminal_nodes) {

    # We can now form the matrices M1 and M2
    curr_M1 <- M1[if(j==1) which_left else which_right,]
    curr_M2 <- M2[if(j==1) which_left else which_right,]
    curr_df <- df[if(j==1) which_left else which_right,]
    curr_n <- nrow(curr_df)

    # We need these matrices:
    n_phi <- nlevels(curr_df$continent)
    n_gamma <- nlevels(curr_df$country)
    M1M1T <- tcrossprod(curr_M1)
    M2M2T <- tcrossprod(curr_M2)
    Delta_R <- diag(curr_n)/tau + M1M1T/(P * tau_gamma)
    Phi_R <- Delta_R + M2M2T/(P * tau_phi)

    # Need the precision matrices
    prec_mu <- sum(solve(Phi_R)) + tau_mu
    prec_phi <- t(curr_M2)%*%solve(Delta_R, curr_M2) + (P*tau_phi) * diag(n_phi)
    prec_gamma <- tau * crossprod(curr_M1) + tau_gamma * P * diag(n_gamma)

    # Update mu
    mean <- sum(solve(Phi_R, curr_df$R)) / prec_mu
    mu_store[[j]][i] <- rnorm(1, mean = mean, sd = sqrt( 1 / prec_mu ))

    # Update phi
    mean <- solve(prec_phi, (t(curr_M2)%*%solve(Delta_R, curr_df$R) + (P*tau_phi * mu_store[[j]][i])))
    # phi_store[[j]][i,] <- mvnfast::rmvn(1, mu = mean,
    #                                     sigma = solve(prec_phi))
    phi_store[[j]][i,] <- rnorm(n_groups, mean = mean,
                                sd = sqrt(1/diag(prec_phi)))

    # Update gamma
    mean <- solve(prec_gamma, tau * t(curr_M1) %*% curr_df$R + tau_gamma * P * phi_store[[j]][i,nesting])
    # gamma_store[[j]][i,] <- mvnfast::rmvn(1, mu = mean,
    #                                       sigma = solve(prec_gamma))
    gamma_store[[j]][i,] <- rnorm(length(nesting), mean = mean,
                                  sd = sqrt(1/diag(prec_gamma)))
  }
}
stop()
# Create some plots of the output -----------------------------------------

# Mu
j <- 1 # Choose which terminal node
data.frame(mu = mu_store[[j]]) |>
  ggplot(aes(x = mu)) +
  geom_histogram() +
  xlim(-3,3) +
  geom_vline(xintercept = true_mu[j], col = "red")

# Phi - re-run these to get different plots
pick <- sample(1:length(true_phi[[j]]), 1)
data.frame(phi = phi_store[[j]][,pick]) |>
  ggplot(aes(x = phi)) +
  geom_histogram() +
  labs(title = pick) +
  xlim(-3,3) +
  geom_vline(xintercept = true_phi[[j]][pick],
             col = 'red')

# Gamma - re-run these to get different plots
pick <- sample(1:length(nesting), 1)
data.frame(gamma = gamma_store[[j]][,pick]) |>
  ggplot(aes(x = gamma)) +
  geom_histogram() +
  labs(title = pick) +
  xlim(-4,4) +
  geom_vline(xintercept = unlist(true_gamma[[j]])[pick],
             col = 'red')
#-----------------------------------------------------------------------
hb_model <- nhebart(formula,
                    data           = df,
                    group_variable = group_variable,
                    subgroup_variable = subgroup_variable,
                    num_trees = 2,
                    priors = list(
                      alpha = 0.95, # Prior control list
                      beta = 2,
                      nu = 2,
                      lambda = 0.1,
                      tau_mu = 16 * num_trees,
                      shape_sigma_phi = 0.5,
                      scale_sigma_phi = 1,
                      sample_sigma_phi = TRUE,
                      shape_sigma_gamma = 0.5,
                      scale_sigma_gamma = 1
                    ), 
                    inits = list(tau         = 1,
                                 tau_phi     = 3, 
                                 sigma_phi   = 1,
                                 sigma_gamma = 1),
                    MCMC = list(iter = 250, # Number of iterations
                                burn = 100, # Size of burn in
                                thin = 1,
                                sigma_phi_sd   = 2,
                                sigma_gamma_sd = 2)
)
#hb_model$trees[[10]]
hb_model
hb_model$rmse
hb_model$mse
#1/hb_model$tau
1/sqrt(hb_model$sigma_gamma)
1/sqrt(hb_model$sigma_phi)

#1/sqrt(hb_model$sigma_gamma)
#1/sqrt(hb_model$sigma_gamma)
pp <- predict_nhebart(newX = df, 
                      new_groups = df$continent, 
                      new_subgroups = df$country,
                      hebart_posterior = hb_model, 
                      type = "mean")
sqrt(mean((pp - df$R)^2)) # 1.907164
cor(pp, df$R)  #


pick <- 2
mus <- vector(length = length(hb_model$trees))
for(k in 2:length(hb_model$trees)){
  mus[k] <- hb_model$trees[[k]][[1]]$tree_matrix[, "mu"][2]
}
mean(mus)
mean(mu_store[[1]])
# MUS -> SAME

pick <- 2
gammas <- vector(length = length(hb_model$trees))
for(k in 2:length(hb_model$trees)){
  gammas[k] <- hb_model$trees[[k]][[1]]$tree_matrix[, "gamma1_3"][3]
}
mean(gammas)
unlist(true_gamma[[2]])[3]


phis <- vector(length = length(hb_model$trees))
for(k in 2:length(hb_model$trees)){
  phis[k] <- hb_model$trees[[k]][[1]]$tree_matrix[, "phi2"][3]
}
mean(phis)
unlist(true_phi[[2]])


data.frame(phi = phi_store[[j]][,pick]) |>
  ggplot(aes(x = phi)) +
  geom_histogram() +
  labs(title = pick) +
  xlim(-3,3) +
  geom_vline(xintercept = true_phi[[j]][pick],
             col = 'red')


data.frame(gamma = gamma_store[[j]][,pick]) |>
  ggplot(aes(x = gamma)) +
  geom_histogram() +
  labs(title = pick) +
  xlim(-4,4) +
  geom_vline(xintercept = unlist(true_gamma[[j]])[pick],
             col = 'red')
