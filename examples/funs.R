get_group_predictions <- function(trees, X, groups, 
                                  subgroups = NULL, 
                                  single_tree = FALSE, 
                                  old_groups) {
  
  
  train_groups <- unique(old_groups)
  new_groups <- unique(groups)
  # Are those new groups?
  which_new <- new_groups[!(new_groups %in% train_groups)]
  if(length(which_new) > 1){
    X_new <- X[groups %in% which_new, ]
    inds_new <- which(groups %in% which_new)
    inds_old <- which(!groups %in% which_new)
    pred_new <- get_predictions(trees, X_new, single_tree = single_tree)
    X_old <- X[!(groups %in% which_new), ]
    if(nrow(X_new) == nrow(X)){
      return(pred_new)
    }
  } else{
    X_old <- X
  }
  
  # Stop nesting problems in case of multiple trees
  if (is.null(names(trees)) & (length(trees) == 1)) trees <- trees[[1]]
  
  if(is.null(subgroups)){
    group_names     <- unique(groups)
    num_groups      <- length(unique(groups))
    group_col_names <- paste0("phi", group_names)
    group_col_names_all <- paste0("phi", groups)
    
    # Normally trees will be a list of lists but just in case
    if (single_tree) {
      # Deal with just a single tree
      if (nrow(trees$tree_matrix) == 1) {
        predictions <- trees$tree_matrix[1, group_col_names][group_col_names_all]
      } else {
        # Loop through the node indices to get predictions
        predictions <- rep(NA, nrow(X_old))
        unique_node_indices <- unique(trees$node_indices)
        # Get the node indices for the current X matrix
        curr_X_node_indices <- fill_tree_details(trees, X_old)$node_indices
        actual_node_indices  <- unique(curr_X_node_indices)
        # Now loop through all node indices to fill in details
        for (i in 1:length(actual_node_indices)) {
          curr_groups <- groups[curr_X_node_indices == actual_node_indices[i]]
          predictions[curr_X_node_indices == actual_node_indices[i]] <-
            trees$tree_matrix[actual_node_indices[i], paste0("phi", curr_groups)]
        }
      }
      # More here to deal with more complicated trees - i.e. multiple trees
    } else {
      # Do a recursive call to the function
      partial_trees <- trees
      partial_trees[[1]] <- NULL # Blank out that element of the list
      predictions <- get_group_predictions(trees[[1]], X_old, groups, single_tree = TRUE, old_groups = old_groups) +
        get_group_predictions(partial_trees, X_old, groups,
                              single_tree = length(partial_trees) == 1,
                              old_groups = old_groups
        )
    }
    
    if(exists("pred_new")){
      final <- data.frame(ind = c(inds_new, inds_old), 
                          predictions = pred_new, predictions)
      predictions <- final$predictions
    }
  } else {
    # With subgroups
    group_names     <- unique(groups)
    num_groups      <- length(unique(groups))
    group_col_names <- paste0("phi", group_names)
    group_col_names_all <- paste0("phi", groups)
    
    subgroup_names     <- unique(subgroups)
    num_subgroups      <- length(unique(subgroups))
    subgroup_col_names <- paste0("gamma", subgroup_names)
    subgroup_col_names_all <- paste0("gamma", subgroups)
    
    
    # Normally trees will be a list of lists but just in case
    if (single_tree) {
      # Deal with just a single tree
      if (nrow(trees$tree_matrix) == 1) {
        vars <- subgroup_col_names
        predictions <- trees$tree_matrix[1, vars][subgroup_col_names_all]
      } else {
        # Loop through the node indices to get predictions
        predictions <- rep(NA, nrow(X_old))
        unique_node_indices <- unique(trees$node_indices)
        # Get the node indices for the current X matrix
        curr_X_node_indices <- fill_tree_details(trees, X_old)$node_indices
        actual_node_indices  <- unique(curr_X_node_indices)
        # Now loop through all node indices to fill in details
        for (i in 1:length(actual_node_indices)) {
          curr_subgroups <- subgroups[curr_X_node_indices == actual_node_indices[i]]
          predictions[curr_X_node_indices == actual_node_indices[i]] <-
            trees$tree_matrix[actual_node_indices[i], paste0("gamma", curr_subgroups)]
        }
      }
      
    } else {
      # Do a recursive call to the function
      partial_trees <- trees
      partial_trees[[1]] <- NULL # Blank out that element of the list
      predictions <- get_group_predictions(trees = trees[[1]], X = X_old, 
                                           groups, subgroups, single_tree = TRUE, old_groups = old_groups) +
        get_group_predictions(partial_trees, X_old, groups, subgroups, 
                              single_tree = length(partial_trees) == 1,
                              old_groups = old_groups
        )
    }
    
    if(exists("pred_new")){
      final <- data.frame(ind = c(inds_new, inds_old), 
                          predictions = pred_new, predictions)
      predictions <- final$predictions
    }
  }
  
  return(predictions)
}


create_stump <- function(num_trees,
                         groups,
                         subgroups, 
                         y,
                         X) {
  
  # Each tree has 8+num_groups columns and 2 elements
  # The 3 elements are the tree matrix, and the node indices
  # The tree matrix has columns:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu values
  # phi values for each group
  # Node size
  
  group_names     <- unique(groups)
  num_groups      <- length(unique(groups))
  
  subgroup_names         <- unique(subgroups)
  #all_subgroup_names     <- paste0("phi_", unique(subgroups))
  num_subgroups          <- length(unique(subgroups))
  
  # Create holder for trees
  all_trees <- vector("list", length = num_trees)
  # Loop through trees
  for (j in 1:num_trees) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] <- vector("list", length = 2)
    # Give the elements names
    names(all_trees[[j]]) <- c(
      "tree_matrix",
      "node_indices"
    )
    # Create the two elements: first is a matrix
    all_trees[[j]][[1]] <- matrix(NA, ncol = 8 + num_groups + num_subgroups, nrow = 1)
    
    # Second is the assignment to node indices
    all_trees[[j]][[2]] <- rep(1, length(y))
    
    # Create column names
    colnames(all_trees[[j]][[1]]) <- c(
      "terminal",
      "child_left",
      "child_right",
      "parent",
      "split_variable",
      "split_value",
      "mu",
      paste0("phi",   sort(group_names)),
      paste0("gamma", sort(subgroup_names)),
      "node_size"
    )
    
    # Set values for stump
    all_trees[[j]][[1]][1, ] <- c(1, 1, NA, NA, NA, NA, 
                                  rep(0, num_groups + num_subgroups + 1), length(y))
  } # End of loop through trees
  
  return(all_trees)
} # End of function


grow_tree <- function(X, y, num_groups, curr_tree, node_min_size) {
  
  
  # Set up holder for new tree
  new_tree <- curr_tree
  
  # Get the list of terminal nodes
  terminal_nodes <- which(new_tree$tree_matrix[, "terminal"] == 1) # Create the list of terminal nodes
  
  # Find terminal node sizes
  terminal_node_size <- new_tree$tree_matrix[terminal_nodes, "node_size"]
  
  # Add two extra rows to the tree in question
  new_tree$tree_matrix <- rbind(
    new_tree$tree_matrix,
    c(1, NA, NA, NA, NA, NA, rep(NA, num_groups + 1), NA), # Make sure they're both terminal
    c(1, NA, NA, NA, NA, NA, rep(NA, num_groups + 1), NA)
  )
  
  # Choose a random terminal node to split
  if(length(terminal_nodes) == 1){
    node_to_split <- terminal_nodes
  } else {
    node_to_split <- sample(stats::na.omit(terminal_nodes), 1, 
                            prob = as.integer(terminal_node_size > node_min_size)
    ) 
    
    ##) # Choose which node to split, set prob to zero for any nodes that are too small
  }
  # Choose a split variable uniformly from all columns
  
  split_value <- 0.5
  split_variable <- 2
  curr_parent <- new_tree$tree_matrix[node_to_split, "parent"] # Make sure to keep the current parent in there. Will be NA if at the root node
  new_tree$tree_matrix[node_to_split, 1:6] <- c(
    0, # Now not temrinal
    nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
    nrow(new_tree$tree_matrix), # child_right is penultimate row
    curr_parent,
    split_variable,
    split_value
  )
  
  #  Fill in the parents of these two nodes
  new_tree$tree_matrix[nrow(new_tree$tree_matrix), "parent"]     <- node_to_split
  new_tree$tree_matrix[nrow(new_tree$tree_matrix) - 1, "parent"] <- node_to_split
  
  # Now call the fill function on this tree
  new_tree <- fill_tree_details(new_tree, X)
  
  # Reject this tree is any values are smaller the node_min_size
  if (any(new_tree$tree_matrix[, "node_size"] < node_min_size)) new_tree <- curr_tree
  
  # Return new_tree
  return(new_tree)
} # End of grow_tree function



simulate_mu_hebart <- function(tree, R, 
                               tau, tau_mu, tau_phi, tau_gamma, 
                               M_1, M_2, num_trees) {
  
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal     <- which(tree$tree_matrix[, "terminal"] == 1)
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  for (i in 1:length(nj)) {
    M_1_j <- M_1[tree$node_indices == which_terminal[i], , drop = FALSE]
    M_2_j <- M_2[tree$node_indices == which_terminal[i], , drop = FALSE]
    R_j <- R[tree$node_indices == which_terminal[i], drop = FALSE]
    
    # Psi_R OK
    #Psi_R <- (diag(nj[i])/tau + 
    #            tcrossprod(M_1_j)/(num_trees*tau_gamma) +
    #            tcrossprod(M_2_j)/(num_trees*tau_phi))
    
    # ANDREW, 23/03/2023
    M1M1T <- tcrossprod(M_1_j)
    M2M2T <- tcrossprod(M_2_j)
    Delta_R <- diag(nj[i])/tau + M1M1T/(num_trees * tau_gamma)
    Psi_R <- Delta_R + M2M2T/(num_trees * tau_phi)
    
    prec <- sum(solve(Psi_R)) + tau_mu
    mean <- sum(solve(Psi_R, R_j)) / prec
    
    # INSTEAD OF:
    # ones <- rep(1, nj[i])
    # Prec_bit <- t(ones)%*%solve(Psi_R, ones) + tau_mu
    # mean <- t(ones)%*%solve(Psi_R, R_j) / Prec_bit
    
    tree$tree_matrix[which_terminal[i], "mu"] <- stats::rnorm(1,
                                                              mean,
                                                              sd = 1/sqrt(prec))
  }
  tree$tree_matrix[which_non_terminal, "mu"] <- NA
  
  return(tree)
}


simulate_phi_hebart <- function(tree, R, groups, 
                                tau, tau_mu, tau_phi, tau_gamma, 
                                M_1, M_2, num_trees) {
  
  # Simulate the group mu values for a given tree
  group_names <- unique(groups)
  num_groups <- length(unique(groups))
  group_col_names <- paste0("phi", group_names)
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Get the group means in each terminal node
  # Doing this with loops but probably can be faster
  for (i in 1:length(nj)) {
    curr_R <- R[tree$node_indices == which_terminal[i]]
    curr_mu <- tree$tree_matrix[which_terminal[i], "mu"]
    curr_M2 <- M_2[tree$node_indices == which_terminal[i], , drop = FALSE]
    curr_M1 <- M_1[tree$node_indices == which_terminal[i], , drop = FALSE]
    delta_R <- diag(nj[i])/tau + tcrossprod(curr_M1)/(num_trees*tau_gamma)
    n_groups <- colSums(curr_M2)
    
    #Prec_one <- solve(t(curr_M2)%*%solve(delta_R, curr_M2) + (num_trees*tau_phi*diag(num_groups)))
    #Prec_two <- (tau*t(curr_M2) %*% curr_R + (num_trees*tau_phi)*curr_mu)
    
    #Prec_one <- solve(t(curr_M2)%*%solve(delta_R, curr_M2) + (num_trees*tau_phi*diag(num_groups)))
    #Prec_two <- (t(curr_M2) %*% solve(delta_R, curr_R) + (num_trees*tau_phi)*curr_mu)
    #mean <- Prec_one %*% Prec_two             
    
    # Same as before
    Prec_one <- t(curr_M2)%*%solve(delta_R, curr_M2) + (num_trees*tau_phi) * diag(num_groups)
    mean <- solve(Prec_one, (t(curr_M2)%*%solve(delta_R, curr_R) + (num_trees*tau_phi*curr_mu)))
    
    
    tree$tree_matrix[which_terminal[i], 
                     sort(group_col_names)] <- mvnfast::rmvn(1,
                                                             mu = mean[, 1],
                                                             sigma = solve(Prec_one))
    
  }
  
  tree$tree_matrix[which_non_terminal, sort(group_col_names)] <- NA
  
  return(tree)
}


simulate_gamma_hebart <- function(tree, R, groups, subgroups, 
                                  tau, tau_mu, tau_phi, tau_gamma, 
                                  M_1, M_2, num_trees) {
  
  
  # Simulate the group mu values for a given tree
  group_names <- unique(groups)
  num_groups <- length(unique(groups))
  group_col_names <- paste0("phi", group_names)
  
  # Simulate the group values for a given tree
  subgroup_names <- unique(subgroups)
  num_subgroups <- length(unique(subgroups))
  subgroup_col_names <- paste0("gamma", subgroup_names)
  
  
  # First find which rows are terminal nodes
  which_terminal <- which(tree$tree_matrix[, "terminal"] == 1)
  which_non_terminal <- which(tree$tree_matrix[, "terminal"] == 0)
  
  # Get node sizes for each terminal node
  nj <- tree$tree_matrix[which_terminal, "node_size"]
  
  # Get the group means in each terminal node
  # Doing this with loops but probably can be faster
  for (i in 1:length(nj)) {
    for(name in group_names){
      which_group <- which(groups == name)
      curr_R <- R[which_group]
      curr_M1 <- M_1[which_group, , drop = FALSE]
      # Get only the columns within that group
      which_M1 <-  colnames(curr_M1) [colSums(curr_M1) > 0]
      curr_M1 <- curr_M1[, colSums(curr_M1) > 0]
      group_name <- paste0("phi", name)
      curr_phi <- tree$tree_matrix[which_terminal[i], group_name]
      
      if(is.vector(curr_M1)){
        curr_M1 <- matrix(curr_M1)
        colnames(curr_M1) <- which_M1
      }
      
      #Prec_one <- solve((tau * t(curr_M1) %*% curr_M1) + (tau_gamma * num_trees * diag(ncol(curr_M1))))
      #Prec_two <- ((tau * t(curr_M1) %*% curr_R) + (tau_gamma * num_trees * curr_phi))
      #mean <- Prec_one %*% Prec_two    
      
      prec <- tau * crossprod(curr_M1) + tau_gamma * num_trees * diag(ncol(curr_M1))
      mean <- solve(prec, tau * t(curr_M1) %*% curr_R + tau_gamma * num_trees * curr_phi)
      
      subgroup_col_names_n <- stringr::str_remove_all(colnames(curr_M1), "factor\\(subgroups\\)")
      subgroup_col_names_n <- paste0("gamma", subgroup_col_names_n)
      
      tree$tree_matrix[which_terminal[i], 
                       sort(subgroup_col_names_n)] <- mvnfast::rmvn(1,
                                                                    mu = mean[, 1],
                                                                    sigma = solve(prec))
      
    } 
  }
  
  tree$tree_matrix[which_non_terminal, sort(subgroup_col_names)] <- NA
  
  return(tree)
}

