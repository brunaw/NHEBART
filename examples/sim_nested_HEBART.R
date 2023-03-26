library(tidyverse)
devtools::load_all(".")
# Simulate some data from a nested HEBART model

# I'm going to assume the data are simulated from a sum of two trees
# with 4 main groups and (3,2,5,3) nested groups inside each main group
# The first tree has 3 terminal nodes and splits on x1 < 0
# at the first split and then x2 < 0 on those observations
# on the left hand side
# The second tree just has one split at x2 < 0.5

set.seed(123)
n_obs <- 500
n_groups <- 4
n_groups_nested <- c(3, 2, 5, 3)
n_trees <- 2
n_terminal_nodes <- c(3, 2)
x1 <- runif(n_obs)
x2 <- runif(n_obs)
group <- sample(1:n_groups, size = n_obs, replace = TRUE)
sub_group <- rep(NA, length = length(group))
for(i in 1:length(group)) {
  sub_group[i] <- sample(1:n_groups_nested[group[i]], 
                         size = 1)
}

# Might want to fiddle around with these values a bit
tau_mu <- 3
tau_phi <- 3 
tau_lambda <- 2
tau <- 1

# Simulate terminal node parameters
mu <- vector('list', length = n_trees)
for (i in 1:length(mu)) {
  mu[[i]] <- rnorm(n_terminal_nodes[i], 0, 1/sqrt(tau_mu))  
}

#[2]]
#[1] -0.2990863  0.5357511

# Then for each terminal node there are n_groups phi parameters
# And also each sub-terminal node as well!
phi <- lambdas <- vector('list', length = n_trees)
for (i in 1:length(phi)) {
  phi[[i]] <- vector('list', length = n_terminal_nodes[i])
  lambdas[[i]] <- vector('list', length = n_terminal_nodes[i])
  for(j in 1:length(phi[[i]])) {
    
    phi[[i]][[j]] <- rnorm(n_groups, mu[[i]][j], sqrt(n_trees / tau_phi))
    lambdas[[i]][[j]] <- vector('list', length = n_groups)
    for (k in 1:length(phi[[i]][[j]])) {
      lambdas[[i]][[j]][[k]] <- rnorm(n_groups_nested[k], 
                                  phi[[i]][[j]][k], 
                                  sqrt(n_trees / tau_lambda))
    }
  }
}
# For phi the order of the indices is tree, terminal node, group
# So phi[[2]][[1]][3] is the 2nd tree, 1st terminal node, 3rd group
# For lambda it's the same but with sub-group added on 
# lambda[[1]][[2]][[2]][3] is the 1st tree, 2nd terminal node, 2nd group, 
# 3rd sub-group
# wow!

# Write the functions for each of the trees
tree_fun1 <- function(x1, x2, group, subgroup) {
  ans <- if(x1 < 0 & x2 < 0) {
    lambda[[1]][[1]][[group]][subgroup]
  } else if(x1 < 0 & x2 >= 0) {
    lambdas[[1]][[2]][[group]][subgroup]
  } else {
    lambdas[[1]][[2]][[group]][subgroup]
  }
  return(ans)
}

tree_fun2 <- function(x1, x2, group, subgroup) {
  ans <- if(x2 < 0.5 ) {
    lambdas[[2]][[1]][[group]][subgroup]
  } else {
    lambdas[[2]][[2]][[group]][subgroup]
  }
  return(ans)
}

# Finally simulate some y values from the tree functions
y <- rep(NA, length = n_obs)

for (i in 1:n_obs) {
  y[i] <- rnorm(1, mean = 
                  tree_fun1(x1[i], x2[i], group[i], sub_group[i]) + 
                  tree_fun2(x1[i], x2[i], group[i], sub_group[i]), 
                sd = 1 / sqrt(tau))
}
# single tree using x2 < 0.5
# Not a particularly helpful plot but just to check something got simulated
library(ggplot2)
df <- data.frame(
  x1 = x1,
  x2 = x2,
  group = group,
  sub_group = sub_group,
  y = y
)
ggplot(df, aes(x = x2, y = y)) +
  geom_point() + 
  facet_grid(group ~ sub_group)


data <- data.frame(
  y = y, 
  x1 = x1, 
  x2 = x2, 
  group = group, 
  subgroup = sub_group
) |> 
  mutate(subgroup = paste0(group, subgroup))


group_variable = "group"
subgroup_variable = "subgroup"
formula = y ~ x1 + x2
num_trees <- 6

tau         = 1
tau_phi     = 3 
tau_mu      = 3
tau_gamma   = 2

hb_model <- nhebart(formula,
                    data           = data,
                    group_variable = group_variable,
                    subgroup_variable = subgroup_variable,
                    num_trees = num_trees,
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
                    MCMC = list(iter =200, # Number of iterations
                                burn = 50, # Size of burn in
                                thin = 1,
                                sigma_phi_sd   = 2,
                                sigma_gamma_sd = 2)
)
#hb_model$trees[[10]]
hb_model
hb_model$rmse
hb_model$mse
1/hb_model$tau

1/sqrt(hb_model$sigma_gamma)
#1/sqrt(hb_model$sigma_gamma)
pp <- predict_nhebart(newX = data, 
                      new_groups = data$group, 
                      new_subgroups = data$subgroup,
                      hebart_posterior = hb_model, 
                      type = "mean")
sqrt(mean((pp - data$y)^2)) # 1.907164
cor(pp, data$y)  #
#----------------------------------------------------------------------  
lme_ss <- lme4::lmer(y ~ x1 + x2 + (1|group) + (1|group:subgroup), data)
pplme <- predict(lme_ss, data)
sqrt(mean((pplme - data$y)^2)) # 3.991818
cor(pplme, data$y) # 0.936175

  
  
