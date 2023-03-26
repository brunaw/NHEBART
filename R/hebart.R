#' @author Bruna Wundervald, \email{brunadaviesw@gmail.com}, Andrew Parnell
#' @export
#' @title Hierarchical Embedded Bayesian Additive Regression Trees
#' @description This function runs a hebart model and returns the tree and
#' other results obtained in the last iteration of the MCMC
#' @param formula The model formula
#' @param data The data to be used in the modeling
#' @param group_variable The name of the grouping variable
#' @param num_trees The number of trees (P)
#' @param control A list with control settings
#' @param priors A list with prior hyperparameters as defined by the model
#' @param inits A list with initial values for parameters
#' @param MCMC A list with MCMC parameters
#' @return A list containing:
#'  Everything
#' @details
#' Priors used ----------------------------------------------------------
#' y_{ij} ~ Normal(m_j, tau^-1)
#' tau    ~ Gamma(nu/2, nu*lambda/2)
#' mu     ~ Normal(0, tau_mu^-1)
#' phi    ~ Normal(mu, sigma_phi^2 / T)
#' ----------------------------------------------------------------------

nhebart <- function(formula,
                   data,
                   group_variable,
                   subgroup_variable, 
                   # X is the feature matrix, y is the target,
                   # groups, # groups is the group number of each obs
                   num_trees = 4, # Number of trees
                   control = list(node_min_size = 2), # Size of smallest nodes
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
                   inits = list(
                     tau = 1,
                     tau_phi = 1, 
                     sigma_phi = 1,
                     sigma_gamma = 1
                   ), # Initial values list
                   MCMC = list(
                     iter = 50, # Number of iterations
                     burn = 10, # Size of burn in
                     thin = 1,
                     sigma_phi_sd   = 2,
                     sigma_gamma_sd = 2
                   )) {

  # Handling formula interface
  formula_int <- stats::as.formula(paste(c(formula), "- 1"))
  response_name <- all.vars(formula_int)[1]
  names_x <- all.vars(formula_int[[3]])

  # Used in create_S to avoid error with stumps
  mod.mat <<- function(f) {
    if(nlevels(f) == 1) {
      m <- matrix(1, nrow = length(f), ncol = 1)
    } else {
      m <- stats::model.matrix(~ f - 1)  
    }
    return(m)
  }
  
  #-------------------------------------------------------
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
  nesting <- data |> 
    dplyr::group_by(group, subgroup) |> 
    dplyr::slice(1) |> 
    dplyr::pull(group)
  #-------------------------------------------------------
  # Extract control parameters
  node_min_size <- control$node_min_size

  # Extract hyper-parameters
  alpha  <- priors$alpha   # Tree shape parameter 1
  beta   <- priors$beta     # Tree shape parameter 2
  tau_mu <- priors$tau_mu # Overall mean precision
  shape_sigma_phi  <- priors$shape_sigma_phi # Weibull prior parameters
  scale_sigma_phi  <- priors$scale_sigma_phi
  sample_sigma_phi <- priors$sample_sigma_phi
  
  shape_sigma_gamma  <- priors$shape_sigma_gamma # Weibull prior parameters
  scale_sigma_gamma  <- priors$scale_sigma_gamma

  # Extract initial values
  tau         <- inits$tau
  sigma       <- 1 / sqrt(tau)
  sigma_phi   <- inits$sigma_phi
  sigma_gamma <- inits$sigma_gamma
  tau_phi     <- 1 / (sigma_phi^2)
  tau_gamma   <- 1 / (sigma_gamma^2)
  log_lik     <- 0

  # Extract MCMC details
  iter <- MCMC$iter # Number of iterations
  burn <- MCMC$burn # Size of burn in
  thin <- MCMC$thin # Amount of thinning
  sigma_phi_sd   <- MCMC$sigma_phi_sd # SD parameter for sigma_phi MH update
  sigma_gamma_sd <- MCMC$sigma_gamma_sd # SD parameter for sigma_phi MH update
  
  # Storage containers
  store_size        <- (iter - burn) / thin
  tree_store        <- vector("list", store_size)
  sigma_store       <- rep(NA, store_size)
  y_hat_store       <- matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store     <- rep(NA, store_size)
  sigma_phi_store   <- rep(NA, store_size)
  sigma_gamma_store <- rep(NA, store_size)
  tau_store         <- rep(NA, store_size)
  mse_store         <- rep(NA, store_size)

  # Scale the response target variable
  y_min <- min(y)
  y_max <- max(y)
  y_scale <- (y - y_min)/(y_max - y_min) - 0.5
  n <- length(y_scale)

  # --------------------------------------------------------------------
  # Finding a value for the parameters of the prior to tau
  #---------------------------------------------------------------------
  lme_form <- paste(paste(c(formula)), 
                    "+ (1|", "group", ")",
                    "+ (1|", "group:subgroup", ")")
  lme_form <- stats::as.formula(lme_form)
  data_lme <- dplyr::mutate(data, y = y_scale)
  my_lme <- lmerTest::lmer(lme_form, data_lme)
  #my_lme <- nlme::lme(lme_form, data_lme, random = (~1|group/subgroup))
  res      <- lme4:::sigma.merMod(my_lme)
  
  nu     <- priors$nu         # Parameter 1 for precision
  lambda <- priors$lambda # Parameter 2 for precision
  p_inv  <- invgamma::pinvgamma(q = res, shape = nu/2, rate = nu*lambda/2)
  
  # Putting high probabilities of the BCART improving a linear model
  while(p_inv < 0.95){
    p_inv <- invgamma::pinvgamma(q = res, shape = nu/2, rate = nu*lambda/2)
    if(p_inv < 0.95){
      nu = abs(nu + stats::rnorm(1))
      lambda = abs(lambda + stats::rnorm(1))
    }
  }
  
  # --------------------------------------------------------------------
  # Finding a value for the parameters of the prior to sigma_phi
  #---------------------------------------------------------------------
  random_effect     <- sqrt(as.data.frame(lme4::VarCorr(my_lme))$vcov[2])
  pr <- parameters::model_parameters(my_lme, effects = "random",
                                     ci_random = TRUE,
                                     verbose = FALSE)
  se <- pr$SE[2]
  random_effect_var <- se^2
  # Sometimes the SE is not available? 
  
  if(!is.na(se) & se > 0){
    shape_sigma_phi  <-  (random_effect^2)/random_effect_var
    scale_sigma_phi  <-  random_effect/random_effect_var
  }
  # Correct parametrization ---------
  scale_sigma_phi <- 1/scale_sigma_phi

  # --------------------------------------------------------------------
  # Finding a value for the parameters of the prior to sigma_gamma
  #---------------------------------------------------------------------
  random_effect     <- sqrt(as.data.frame(lme4::VarCorr(my_lme))$vcov[1])
  se <- pr$SE[1] # phi:gamma
  random_effect_var <- se^2
  
  if(!is.na(se) & se > 0){
    shape_sigma_gamma  <-  (random_effect^2)/random_effect_var
    scale_sigma_gamma  <-  random_effect_var/random_effect
  }
  # Correct parametrization ---------
  scale_sigma_gamma <- 1/scale_sigma_gamma
  
  #---------------------------------------------------------------------
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
    groups = groups,
    subgroups = subgroups, 
    y = y_scale,
    X = X
  )
  predictions <- get_group_predictions(
    trees  = curr_trees, X, groups, subgroups, 
    single_tree = num_trees == 1, 
    old_groups = groups)
  
  # Set up a progress bar
  pb <- utils::txtProgressBar(
    min = 1, max = iter,
    style = 3, width = 60,
    title = "Running rBART..."
  )
  
  # Start the iterations loop
  for (i in 1:iter) {
    utils::setTxtProgressBar(pb, i)

    
    # If at the right place store everything
    if ((i > burn) & ((i %% thin) == 0)) {
      curr <- (i - burn) / thin
      tree_store[[curr]]  <- curr_trees
      sigma_store[curr]   <- sigma
      y_hat_store[curr, ] <- predictions
      log_lik_store[curr] <- log_lik
      sigma_phi_store[curr]   <- sigma_phi
      sigma_gamma_store[curr] <- sigma_gamma
      tau_store[curr] <- tau
      mse_store[curr] <- mse 
    }
 
    
    # Start looping through trees
    for (j in 1:num_trees) {

      # Calculate partial residuals for current tree
      if (num_trees > 1) {
        partial_trees <- curr_trees
        partial_trees[[j]] <- NULL # Blank out that element of the list
        current_partial_residuals <- y_scale -
          get_group_predictions(
            trees = partial_trees, X, groups, subgroups, 
            single_tree = num_trees == 2, 
            old_groups = groups
          )
      } else {
        current_partial_residuals <- y_scale
        #current_partial_residuals <- y
      }
      
      
      # Propose a new tree via grow/change/prune/swap
      #type <- sample(c("grow", "prune", "change", "swap"), 1)
      if (i < max(floor(0.1 * burn), 10)) type <- "grow" # Grow for the first few iterations

      type <- sample(c("grow", "prune"), 1, prob = c(0.50, 0.50))
      # Get a new tree!
      new_trees <- curr_trees

      new_trees[[j]] <- update_tree(
        y = y_scale,
        X = X,
        groups = groups,
        subgroups = subgroups,
        type = type,
        curr_tree = curr_trees[[j]],
        node_min_size = node_min_size
      )
      
      # Calculate the complete conditional and acceptance probability
      l_new <- full_conditional_nhebart(
        tree = new_trees[[j]],
        R = current_partial_residuals,
        num_trees,
        tau,
        tau_mu, 
        tau_phi,
        tau_gamma,
        M_1, 
        M_2
      ) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      l_old <- full_conditional_nhebart(
        curr_trees[[j]],
        current_partial_residuals,
        num_trees,
        tau,
        tau_mu, 
        tau_phi,
        tau_gamma,
        M_1, 
        M_2
      ) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # If accepting a new tree update all relevant parts
      log.alpha <- (l_new - l_old)
      accept <- log.alpha >= 0 || log.alpha >= log(stats::runif(1))
      if (accept) curr_trees <- new_trees
      
      # Update mu whether tree accepted or not
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
        num_trees, 
        nesting = nesting
      )
      

      # Check the trees
      if (any(curr_trees$tree_matrix[, "node_size"] < node_min_size)) browser()
    } # End loop through trees

    # Calculate full set of predictions
    predictions <- get_group_predictions(trees = curr_trees, 
                                         X, 
                                         groups, 
                                         subgroups = subgroups, 
                                         single_tree = num_trees == 1,
                                         old_groups = groups
    )

    mse <- mean((y_scale - predictions)^2)
    
    # Update tau
    tau <- update_tau(
      y = y_scale,
      predictions,
      nu,
      lambda
    )
    sigma <- 1 / sqrt(tau)
    
    S1 <- create_S(curr_trees, groups = subgroups)
    S2 <- create_S(curr_trees, groups = groups)
    
    # Update tau_phi
    if(sample_sigma_phi){
      sigma_phi <- update_sigma_phi(
        y_scale, sigma_phi, tau_mu, tau, tau_gamma, 
        shape_sigma_phi, scale_sigma_phi, 
        num_trees, sigma_phi_sd, S1, S2
      )
    }
    
    tau_phi <- 1 / (sigma_phi^2)
    
    sigma_gamma <- update_sigma_gamma(
      y_scale, sigma_gamma, tau_phi, tau_mu, tau, 
      shape_sigma_gamma, scale_sigma_gamma, num_trees, 
      sigma_gamma_sd = 0.3, S1, S2
    )
    tau_gamma <- 1 / (sigma_gamma^2)
  
  
    # The end ---------

    # Get the overall log likelihood -- UPDATE HERE
    Omega_y <- calculate_omega(num_trees, tau, tau_phi = tau_phi, tau_gamma = tau_gamma, tau_mu, S1, S2, n)
    log_lik <- mvnfast::dmvn(y, rep(0, n), Omega_y, log = TRUE)
      
  } # End iterations loop
  cat("\n") # Make sure progress bar ends on a new line

  result <- list(
    trees = tree_store,
    sigma = sigma_store,
    y_hat = (y_hat_store + 0.5) * (max(y) - min(y)) + min(y),
    log_lik = log_lik_store,
    sigma_phi   = sigma_phi_store,
    sigma_gamma = sigma_gamma_store,
    tau = tau_store, 
    mse = mse_store, 
    y = y,
    X = X,
    groups = unique(groups), 
    iter = iter,
    burn = burn,
    thin = thin,
    store_size = store_size,
    num_trees = num_trees,
    formula = formula,
    y_min = y_min,
    y_max = y_max
  )

  # RMSE calculation
  pred <- predict_nhebart(newX = data, new_groups = groups, 
                          new_subgroups = subgroups, 
                         hebart_posterior = result, type = "mean")
  mse <- mean((pred - y)^2)
  
  rmse <- sqrt(mse)
  r.squared <- 1 - mse / stats::var(y)

  result$rmse <- rmse
  result$r.squared <- r.squared
  result$num_variables <- length(names_x)

  class(result) <- "hebart"

  return(result = result)
} # End main function