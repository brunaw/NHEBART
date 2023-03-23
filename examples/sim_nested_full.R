#------------------------------------------------------------------------------
library(tidyverse)
source('examples/funs.R')
# Simulate some data from a nested HEBART model
# I'm going to assume the data are simulated from a sum of two trees
# with 4 main groups and (3,2,5,3) nested groups inside each main group
# The tree just has one split at x2 < 0.5

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
                  tree_fun2(x1[i], x2[i], group[i], sub_group[i]), 
                sd = 1 / sqrt(tau))
}
library(ggplot2)

data <- data.frame(
  y = y, 
  x1 = x1, 
  x2 = x2, 
  group = group, 
  subgroup = sub_group
) |> 
  mutate(subgroup = paste0(group, subgroup))

# ggplot(data, aes(x = x2, y = y)) +
#   geom_point() + 
#   facet_grid(group ~ subgroup)

#------------------------------------------------------------------------------
group_variable = "group"
subgroup_variable = "subgroup"
formula = y ~ x1 + x2
num_trees <- 1

tau         = 1
tau_phi     = 3 
tau_mu      = 3
tau_gamma   = 2
#------------------------------------------------------------------------------
# Setting initial variables
data <- dplyr::select(
  data, c(!!response_name, !!names_x, !!group_variable, !!subgroup_variable))
# data    <- dplyr::select(data, c(!!response_name, !!names_x, "group"))

names(data)[names(data) == group_variable]    <- "group"
names(data)[names(data) == subgroup_variable] <- "subgroup"
groups    <- data$group
subgroups <- data$subgroup
mf <- stats::model.frame(formula_int, data = data)
X <- as.matrix(stats::model.matrix(formula_int, mf))
y <- stats::model.extract(mf, "response")


# Get the group matrix M2 (as in the paper)
M_2 <- stats::model.matrix(~ factor(groups) - 1)
group_sizes <- table(groups)
num_groups <- length(group_sizes)

# Get the subgroup matrix M1
M_1 <- stats::model.matrix(~ factor(subgroups) - 1)
subgroup_sizes <- table(subgroups)
num_subgroups <- length(subgroup_sizes)

# Create a list of trees for the initial stump
curr_trees <- create_stump(
  num_trees = num_trees,
  groups    = groups,
  subgroups = subgroups, 
  y = y,
  X = X
)

predictions <- get_group_predictions(
  trees  = curr_trees, X, groups, subgroups, 
  single_tree = num_trees == 1, 
  old_groups = groups)

#------------------------------------------------------------------------------
iters <- 50
tree_store        <- vector("list", iters)
mse_store         <- rep(NA, iters)
corr_store        <- rep(NA, iters)
# Just do one single tree growing at x2 > 0.

current_partial_residuals = y
j = 1 # only one tree anyways

for (i in 1:iter) {

  curr_trees[[j]] <- grow_tree(
    y = y,
    X = X,
    curr_tree = curr_trees[[j]], 
    num_groups =num_subgroups + num_groups, 
    node_min_size = 1
  )
  
 
  curr_trees[[j]] <- simulate_mu_hebart(
    tree = curr_trees[[j]],
    R = current_partial_residuals,
    tau,
    tau_mu, 
    tau_phi,
    tau_gamma,
    M_1, 
    M_2,
    num_trees
  )
  
  # Update phi as well
  curr_trees[[j]] <- simulate_phi_hebart(
    tree = curr_trees[[j]],
    R = current_partial_residuals,
    groups,
    tau,
    tau_mu, 
    tau_phi,
    tau_gamma,
    M_1, 
    M_2,
    num_trees
  )
  
  curr_trees[[j]] <- simulate_gamma_hebart(
    tree = curr_trees[[j]],
    R = current_partial_residuals,
    groups,
    subgroups, 
    tau,
    tau_mu, 
    tau_phi,
    tau_gamma,
    M_1, 
    M_2,
    num_trees
  )
  
  predictions <- get_group_predictions(trees = curr_trees, 
                                       X, 
                                       groups, 
                                       subgroups = subgroups, 
                                       single_tree = num_trees == 1,
                                       old_groups = groups
  )
  
  mse <- mean((y - predictions)^2)
  mse_store[i] <- mse
  corr_store[i] <- cor(y, predictions)
  tree_store[[i]] <- curr_trees
}



