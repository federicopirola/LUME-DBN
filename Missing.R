B_Init <- function(df, init = 'zeros', sample_id = 'Sample_Id', time_id = 'Time') {
  # Extract variable names (excluding sample_id, time_id, intercept column)
  variable_names <- setdiff(colnames(df), c(sample_id, time_id, "1"))
  
  k <- length(variable_names)
  
  # Create named matrix: rows = target variables, columns = intercept + predictors
  if (init == 'zeros'){
    B <- matrix(0, nrow = k, ncol = k + 1)
  }
  else{
    B <- init
  }
  rownames(B) <- variable_names
  colnames(B) <- c("1", variable_names)
  
  B
}

M_Init <- function(df, init = 'random', sample_id = 'Sample_Id', time_id = 'Time', intercept_move = FALSE, max_parents = ncol(df) - 2) {
  # Extract variable names (excluding sample_id, time_id, intercept column)
  variable_names <- setdiff(colnames(df), c(sample_id, time_id, "1"))
  
  k <- length(variable_names)
  
  # Create named matrix: rows = target variables, columns = intercept + predictors
  if (init == 'random'){
    #M <- cbind(1, matrix(rbinom(k^2, 1, 0.5), nrow = k))
    M <- matrix(0, nrow = k, ncol = k + 1)
    for (i in 1:k) {
      possible_parents <- setdiff(1:k, i)  # No self-loops
      num_parents <- sample(0:min(max_parents, k - 1), 1)
      selected_parents <- sample(possible_parents, num_parents)
      M[i, selected_parents + 1] <- 1  # +1 accounts for intercept column
    }
    M[, 1] <- 1
  }
  else if (init == 'empty'){
    M <- cbind(ifelse(intercept_move,0,1), matrix(0, nrow = k, ncol = k))
  }
  else if (init == 'complete'){
    M <- matrix(1, nrow = k, ncol = k + 1)
  }
  else{
    M <- init
  }
  rownames(M) <- variable_names
  colnames(M) <- c("1", variable_names)
  
  M
}

InitialMissingImputation <- function(df, s = 5, sample_id = "Sample_Id", time_id = "Time"){
  df_imputed <- df
  
  # Get variable columns (all except ID columns)
  variable_cols <- setdiff(names(df), c(sample_id, time_id))
  
  # Extract time vector
  time_points <- sort(unique(df[[time_id]]))
  
  for (t in time_points) {
    for (var in variable_cols) {
      # Identify rows at time t where var is missing
      idx_missing <- which(df[[time_id]] == t & is.na(df[[var]]))
      if (length(idx_missing) == 0) next
      
      # 1. Mean at time t
      values_t <- df[df[[time_id]] == t, var, drop = TRUE]
      if (sum(!is.na(values_t)) >= s) {
        #impute_val <- mean(values_t, na.rm = TRUE)
        impute_val <- rnorm(1,mean(values_t, na.rm = TRUE),0.1)
        
      } else {
        # 2. Neighbors [t-1, t+1]
        neighbor_times <- c()
        if (t > min(time_points)) neighbor_times <- c(neighbor_times, t - 1)
        if (t < max(time_points)) neighbor_times <- c(neighbor_times, t + 1)
        
        values_neighbors <- df[df[[time_id]] %in% neighbor_times, var, drop = TRUE]
        if (sum(!is.na(values_neighbors)) >= s) {
          #impute_val <- mean(values_neighbors, na.rm = TRUE)
          impute_val <- rnorm(1,mean(values_neighbors, na.rm = TRUE),0.1)
          
        } else {
          # 3. Fallback: mean over full data
          #impute_val <- mean(df[[var]], na.rm = TRUE)
          impute_val <- rnorm(1,mean(df[[var]], na.rm = TRUE),0.1)
        }
      }
      
      # Impute missing values at time t
      df_imputed[idx_missing, var] <- impute_val
    }
  }
  
  df_imputed
}



MissingFCDComponents <- function(B, M, Sigma) {
  var_names <- rownames(M)   # variable names (excluding intercept)
  intercept_name <- colnames(M)[1]  # assuming first column is intercept, e.g. "(Intercept)"
  
  FCDComponents <- list()
  
  for (i in seq_along(var_names)) {
    target_var <- var_names[i]
    
    # Parents of target_var (including intercept)
    parents_i <- colnames(M)[which(M[i, ] == 1)]  # names of parents
    
    # Children variables where target_var is a parent (note: target_var is in rows, parents in columns)
    # So children are rows where M[, target_var] == 1
    children_vars <- rownames(M)[which(M[, target_var] == 1)]
    
    # For each child, get the other parents (excluding target_var), betas and sigma^2
    child_info <- lapply(children_vars, function(child_var) {
      other_parents <- setdiff(colnames(M)[which(M[child_var, ] == 1)], target_var)
      list(
        child = child_var,
        beta_i = B[child_var, target_var],          # beta for target_var as parent of child_var
        beta_other = B[child_var, other_parents],   # betas for other parents
        sigma2_j = Sigma[child_var],
        other_covariate_names = other_parents       # named vector of other parents
      )
    })
    
    # Save components for this variable
    FCDComponents[[target_var]] <- list(
      sigma2_i = Sigma[target_var],
      parents_i = parents_i,
      beta_parents_i = B[target_var, parents_i],
      model_children = child_info
    )
  }
  
  FCDComponents
}




MissingFCDSample <- function(data, structure_i) {
  # data: matrix with rows = times (t-1, t, t+1), cols = variable names (no intercept col)
  sigma2_i <- structure_i$sigma2_i
  
  sum_precisions <- 1 / sigma2_i
  weighted_sum <- sum(data[1, structure_i$parents_i] * structure_i$beta_parents_i) / sigma2_i
  
  for (child in structure_i$model_children) {
    j <- child$child
    mu_minus_i_j <- sum(data[2, child$other_covariate_names] * child$beta_other)  # time t
    x_j_tp1 <- data[3, j]  # time t+1
    
    beta_i <- child$beta_i
    sigma2_j <- child$sigma2_j
    
    sum_precisions <- sum_precisions + (beta_i^2) / sigma2_j
    weighted_sum <- weighted_sum + (beta_i / sigma2_j) * (x_j_tp1 - mu_minus_i_j)
  }
  sigma2_post <- 1 / sum_precisions
  rnorm(1, mean = sigma2_post * weighted_sum, sd = sqrt(sigma2_post))
}

MissingFCDSample_t0 <- function(data, structure_i, mean_0, sd_0, variable_name) {
  # data: matrix with rows = times (t=0, t=1), cols = variable names
  sigma2_i <- sd_0[variable_name]
  mu_0 <- mean_0[variable_name]
  
  sum_precisions <- 1 / sigma2_i
  weighted_sum <- mu_0 / sigma2_i
  
  for (child in structure_i$model_children) {
    j <- child$child
    mu_minus_i_j <- sum(data[2, child$other_covariate_names] * child$beta_other)  # time 1
    x_j_t1 <- data[2, j]  # time 1
    
    beta_i <- child$beta_i
    sigma2_j <- child$sigma2_j
    
    sum_precisions <- sum_precisions + (beta_i^2) / sigma2_j
    weighted_sum <- weighted_sum + (beta_i / sigma2_j) * (x_j_t1 - mu_minus_i_j)
  }
  sigma2_post <- 1 / sum_precisions
  rnorm(1, mean = sigma2_post * weighted_sum, sd = sqrt(sigma2_post))
}

MissingFCDSample_tT <- function(data, structure_i) {
  sigma2_i <- structure_i$sigma2_i
  sum_precisions <- 1 / sigma2_i
  weighted_sum <- sum(data[1, structure_i$parents_i] * structure_i$beta_parents_i) / sigma2_i
  sigma2_post <- 1 / sum_precisions
  rnorm(1, mean = sigma2_post * weighted_sum, sd = sqrt(sigma2_post))
}

MissingValuesUpdate <- function(data, missing_mask, structure, sample_id = 'Sample_Id', time_id = 'Time', mean_0, sd_0) {
  samples <- unique(data[[sample_id]])
  times <- unique(data[[time_id]])
  T_max <- max(times)
  
  # Get variable names for columns excluding sample_id and time_id
  var_names <- setdiff(colnames(data), c(sample_id, time_id))
  
  for (s in samples) {
    for (t in 0:T_max) {
      time_window <- c(t - 1, t, t + 1)
      valid_times <- time_window[time_window >= 0 & time_window <= T_max]
      
      sub_data <- data[data[[sample_id]] == s & data[[time_id]] %in% valid_times, ]
      if (nrow(sub_data) < length(valid_times)) next
      
      X <- as.matrix(sub_data[, var_names, drop = FALSE])
      X <- X[match(valid_times, sub_data[[time_id]]), , drop = FALSE]
      #return(X)
      mask_row <- which(data[[sample_id]] == s & data[[time_id]] == t)
      missing_row <- missing_mask[mask_row, var_names, drop = FALSE]
      missing_vars <- var_names[as.logical(missing_row)]
      if (length(missing_vars) == 0) next
      
      for (var in missing_vars) {
        if (t == 0) {
          sampled_val <- MissingFCDSample_t0(X[1:2, , drop = FALSE], structure[[var]], mean_0, sd_0, variable_name = var)
        } else if (t == T_max) {
          sampled_val <- MissingFCDSample_tT(X[1:2, , drop = FALSE], structure[[var]])
        } else {
          sampled_val <- MissingFCDSample(X, structure[[var]])
        }
        X[which(valid_times == t), var] <- sampled_val
        data[mask_row, var] <- sampled_val
      }
    }
  }
  
  data
}

PosteriorSamplesMatrices <- function(posterior_samples){
  # Get the full set of covariates across all variables
  all_covariates <- sort(unique(unlist(lapply(posterior_samples, function(x) x$covariates))))
  
  # Initialize matrices and vectors
  vars <- names(posterior_samples)
  B <- matrix(0L, nrow = length(vars), ncol = length(all_covariates),
              dimnames = list(vars, all_covariates))
  M <- matrix(0, nrow = length(vars), ncol = length(all_covariates),
              dimnames = list(vars, all_covariates))
  Sigma <- setNames(rep(0, length(vars)), vars)
  Delta <- setNames(rep(0, length(vars)), vars)
  
  # Iterate through each variable (target)
  for (v in vars) {
    post <- posterior_samples[[v]]
    covars <- post$covariates
    
    last_beta <- tail(post$betas, 1)[[1]][[1]]      # character vector
    last_model <- tail(post$models, 1)[[1]]    # character vector
    last_sigma2 <- tail(post$sigma2s, 1)[[1]]  # numeric
    last_delta2 <- tail(post$delta2s, 1)[[1]]  # numeric
    
    # Match covariate presence
    B[v, covars] <- last_beta
    
    # Assign corresponding model string (by covariate order)
    M[v, covars] <- last_model
    
    # Scalars
    Sigma[v] <- last_sigma2
    Delta[v] <- last_delta2
  }
  
  # Output
  list(B = B, M = M, Sigma = Sigma, Delta = Delta)
}

PosteriorModelMatricesOverTime <- function(posterior_samples) {
  all_covariates <- sort(unique(unlist(lapply(posterior_samples, function(x) x$covariates))))
  vars <- names(posterior_samples)
  n_iter <- length(posterior_samples[[1]]$models)
  
  model_matrices <- vector("list", n_iter)
  
  for (iter in 1:n_iter) {
    A <- matrix(0, nrow = length(vars), ncol = length(all_covariates),
                dimnames = list(vars, all_covariates))
    
    for (v in vars) {
      model_iter <- posterior_samples[[v]]$models[[iter]]
      
      # Check it's not null and has names
      if (!is.null(model_iter) && length(model_iter) > 0 && !is.null(names(model_iter))) {
        covariate_names <- names(model_iter)
        matched_covs <- intersect(covariate_names, all_covariates)
        
        # Assign values to corresponding columns in the matrix row
        A[v, matched_covs] <- model_iter[matched_covs]
      }
    }
    model_matrices[[iter]] <- A
  }
  
  return(model_matrices)
}

# Generate full dataset (no missing values)
generate_sample_data <- function(N, T_max, k) {
  sample_ids <- rep(1:N, each = (T_max + 1))
  time_ids <- rep(0:T_max, times = N)
  variable_data <- matrix(rnorm(N * (T_max + 1) * k), ncol = k)
  colnames(variable_data) <- paste0("X", 1:k)
  data <- data.frame(sample_id = sample_ids, time_id = time_ids, variable_data)
  return(data)
}

# Introduce missingness
create_missing_mask <- function(data, k, missing_prob = 0.1) {
  variable_matrix <- as.matrix(data[, paste0("X", 1:k)])
  mask <- matrix(runif(nrow(variable_matrix) * k) < missing_prob,
                 nrow = nrow(variable_matrix), ncol = k)
  return(mask)
}

# Apply missingness
apply_missingness <- function(data, missing_mask, k) {
  data_missing <- data
  for (j in 1:k) {
    data_missing[missing_mask[, j], paste0("X", j)] <- NA
  }
  return(data_missing)
}

AverageM <- function(transition_network){
  all_nodes <- names(transition_network)
  # Initialize binary matrix
  parent_matrix <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes),
                          dimnames = list(all_nodes, all_nodes))
  # Fill matrix: row = target, col = parent
  for (target in names(transition_network)) {
    parents <- setdiff(transition_network[[target]]$parents, "1")
    parent_matrix[target, parents] <- 1
  }
  parent_matrix
}

PrecisionRecall <- function(real_transition_network, inclusion_prob_matrix){
  # Flatten both matrices to vectors for easy comparison
  true_vec <- as.vector(real_transition_network)
  pred_vec <- as.vector(inclusion_prob_matrix)
  
  # Compute confusion matrix components
  TP <- sum(true_vec == 1 & pred_vec == 1)  # True Positives
  FP <- sum(true_vec == 0 & pred_vec == 1)  # False Positives
  FN <- sum(true_vec == 1 & pred_vec == 0)  # False Negatives
  
  # Compute precision and recall
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  cat("Precision:", ifelse(is.nan(precision), NA, precision), "\n")
  cat("Recall:", ifelse(is.nan(recall), NA, recall), "\n")
  
  list(precision, recall)
}

InclusionProbMatrix <- function(inclusion_probs) {
  # Get all unique covariates, excluding "1"
  all_covariates <- sort(setdiff(unique(unlist(lapply(inclusion_probs, names))), "1"))
  all_targets <- names(inclusion_probs)
  
  # Initialize matrix with 0s
  prob_matrix <- matrix(0, nrow = length(all_targets), ncol = length(all_covariates),
                        dimnames = list(all_targets, all_covariates))
  
  # Fill in the inclusion probabilities
  for (target in all_targets) {
    probs <- inclusion_probs[[target]]
    probs <- probs[setdiff(names(probs), "1")]  # remove "1" if present
    covars <- names(probs)
    prob_matrix[target, covars] <- unlist(probs)
  }
  
  return(prob_matrix)
}

list_to_array <- function(model_list) {
  P <- nrow(model_list[[1]])
  Q <- ncol(model_list[[1]])
  T <- length(model_list)
  arr <- array(NA, dim = c(P, Q, T),
               dimnames = list(rownames(model_list[[1]]),
                               colnames(model_list[[1]]),
                               NULL))
  for (t in 1:T) {
    arr[,,t] <- model_list[[t]]
  }
  return(arr)
}


ComputePSRF <- function(adj_array, threshold = 1.1, every = 100, stop_after = NULL) {
  P <- dim(adj_array)[1]
  Q <- dim(adj_array)[2]
  T <- dim(adj_array)[3]
  M <- dim(adj_array)[4]
  
  # Set maximum iteration for evaluation
  max_iter <- if (!is.null(stop_after)) min(T, stop_after) else T
  selected_iters <- seq(every, max_iter, by = every)
  rates <- numeric(length(selected_iters))
  
  # Precompute observed edges (at least one nonzero over all chains & time)
  observed_edges <- apply(adj_array, c(1,2), function(x) any(!is.na(x) & x != 0))
  
  for (idx in seq_along(selected_iters)) {
    t <- selected_iters[idx]
    converged_count <- 0
    observed_count <- 0
    
    for (i in 1:P) {
      for (j in 1:Q) {
        if (!observed_edges[i, j]) next
        
        # Samples up to iteration t, across all chains
        samples <- matrix(NA, nrow = t, ncol = M)
        for (m in 1:M) {
          samples[, m] <- adj_array[i, j, 1:t, m]
        }
        
        chain_means <- colMeans(samples, na.rm = TRUE)
        within_vars <- apply(samples, 2, function(x) if (all(is.na(x))) NA_real_ else var(x, na.rm = TRUE))
        W <- mean(within_vars, na.rm = TRUE)
        B <- if (M > 1) t * var(chain_means, na.rm = TRUE) else 0
        
        if (is.na(W) || W == 0) {
          R_hat <- 1
        } else {
          V_hat <- ((t - 1) / t) * W + (1 / t) * B
          R_hat <- sqrt(V_hat / W)
        }
        
        observed_count <- observed_count + 1
        if (R_hat < threshold) {
          converged_count <- converged_count + 1
        }
      }
    }
    
    rates[idx] <- if (observed_count > 0) converged_count / observed_count else NA
  }
  
  return(data.frame(iteration = selected_iters, rate = rates))
}



RJMCMC_DBN_Missing <- function(df,
                       sample_id = 'Sample_Id',
                       time_id = 'Time',
                       parameters_moves = 'Gibbs',
                       model_moves = 'MH-Marginal',
                       type_model_prior = 'Uniform',
                       type_segment_prior = 'Uniform_card', # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
                       type_coupling = 'piece-wise', # type of parameter coupling ('piece-wise' or 'global')
                       coupling = FALSE, # if parameter coupling should be applied
                       standardized = FALSE,
                       time_lags = 1,
                       self_loop = FALSE,
                       p_MCMC = 1, # probability of doing an MCMC move
                       p_RJMCMC = 0.1,
                       p_ChangePoint = 0.1, # probability of doing a Change Point Detection move
                       missing_update_period = 10, #how many epochs between missing updates
                       alpha = 0.01,
                       beta = 0.01,
                       alpha_xi = 0.01, #xi2 scale parameter
                       beta_xi = 0.01, #xi2 rate parameter
                       a = 0.01,
                       b = 0.01,
                       v = 0.1,
                       r = 1,
                       lambda = 1,
                       rho = 0.1, #geometric probability parameter for segment distance prior
                       iota = 1, # expected value of number of successes of Poisson segments' granularity prior
                       intercept_move = FALSE,
                       max_parents = ncol(df) - 3 + self_loop,
                       max_changepoints = NaN, #max number of changepoints
                       epochs = 1000,
                       init_sigma2 = 1,
                       init_delta2 = 1,
                       init_B = 'zeros',
                       init_M = 'random',
                       init_xi2 = 1, # initial xi2
                       init_allocation_vector = NaN, # initial segment set
                       segmentation = FALSE, # if data segmentation should be applied
                       segmentation_var = time_id, # variable on which data is segmented
                       n_segments = NaN, # Fixed number of segments (if desired)
                       segmentation_vector = NaN, # Fixed segmentation vector
                       log = TRUE,
                       sigma2_in_betas_prior = FALSE
) {
  # Get variable names (assumes all columns except Sample_Id and Time are variables)
  posterior_samples <- list()
  missing_values = list()
  variables <- sort(setdiff(names(df), c(sample_id, time_id)))
  k <- length(variables)
  # Define covariates as all variables at t-1
  if (!self_loop) {
    covars <- lapply(variables, function(x) {c('1',unlist(sapply(1:time_lags, function(lag) paste0(setdiff(variables, x), "_t-", lag))))})
  }
  else{
    covars <- lapply(variables, function(x) {c('1',unlist(sapply(1:time_lags, function(lag) paste0(variables, "_t-", lag))))})
  }
  var_names <- setdiff(names(df), c(sample_id, time_id))
  #mean_0 <- sapply(var_names, function(x) mean(df[df[[time_id]] == 0, x], na.rm = TRUE))
  mean_0 <- sapply(var_names, function(x) {
    time0_vals <- df[df[[time_id]] == 0, x]
    if (any(!is.na(time0_vals))) {
      mean(time0_vals, na.rm = TRUE)
    } else {
      time1_vals <- df[df[[time_id]] == 1, x]
      if (any(!is.na(time1_vals))) {
        mean(time1_vals, na.rm = TRUE)
      } else {
        mean(df[[x]], na.rm = TRUE)
      }
    }
  })
  names(mean_0) <- var_names  # ensures names are retained
  sd_0 <- sapply(var_names, function(x) {1}) #sd(df[df[[time_id]] == 0, x], na.rm = TRUE)
  names(sd_0) <- var_names
  df <- cbind(df[,c(sample_id, time_id)], '1' = 1, df[, sort(setdiff(names(df), c(sample_id, time_id, "1")))])
  missing_mask <- df %>%
    mutate(across(-c(sample_id, time_id), ~ ifelse(is.na(.), TRUE, FALSE)))
  if (standardized){
    df[,variables] <- StandardizeData(df[,variables])
  }
  #df <- TemporalIterativeImpute(df)
  df <-InitialMissingImputation(df, s = 5) #TemporalMice(df, sample_id = sample_id) #
  # Create lagged dataset for DBN
  M <- M_Init(df, init = init_M, sample_id = sample_id, time_id = time_id, intercept_move = intercept_move, max_parents = max_parents - intercept_move)
  B <- B_Init(df, init = init_B, sample_id = sample_id, time_id = time_id)
  Sigma <- rep(init_sigma2, k)
  Delta <- rep(init_delta2, k)
  for (j in 1:k){
    posterior_samples[[variables[j]]] <- list('betas' = list(),'sigma2s' = list(),'delta2s' = list(),'models' = list(),'allocation_vectors' = list(), 'covariates'= gsub("_t-1", "", covars[[j]]))
  }
  # Fit one model for each variable at time t
  for (epoch in 1:round(epochs/missing_update_period)) {
    lagged_data <- LaggedDatasetGeneration2(df[,!(names(df) %in% c('1'))], lag = time_lags, sample_id = sample_id, time_id = time_id)
    lagged_data <- cbind(lagged_data[, 1:2], '1' = 1, lagged_data[, 3:ncol(lagged_data)])
    for (j in 1:k) {
      if (self_loop){
        b_init <- B[j,]
        m_init <- M[j,]
      }
      else{
        b_init <- B[j,-(j+1)]
        m_init <- M[j,-(j+1)]
      }
      outcome_var <- paste0(variables[j], "_t")
      #cat(outcome_var,'\n\n')
      res <- RJMCMC(
        data = lagged_data[,c(sample_id, time_id, covars[[j]], outcome_var)],
        epochs = missing_update_period,
        v = v,
        type_segment_prior = type_segment_prior, # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
        type_coupling = type_coupling, # type of parameter coupling ('piece-wise' or 'global')
        coupling = coupling, # if parameter coupling should be applied
        p_MCMC = p_MCMC, # probability of doing a parameter move
        p_RJMCMC = p_RJMCMC, # probability of doing an model move
        p_ChangePoint = p_ChangePoint, # probability of doing a Change Point Detection move
        alpha = alpha,
        beta =  beta,
        alpha_xi = alpha_xi, #xi2 scale parameter
        beta_xi = beta_xi, #xi2 rate parameter
        a = a,
        b = b,
        r = r,
        lambda = lambda,
        rho = rho, #geometric probability parameter for segment distance prior
        iota = iota, # expected value of number of successes of Poisson segments' granularity prior
        max_parents = max_parents,
        max_changepoints = max_changepoints, #max number of changepoints
        intercept_move = intercept_move,
        outcome = outcome_var,
        covariates = covars[[j]],
        init_model = m_init,
        init_betas = b_init,
        betas_mean = rep(0, (length(covars[[j]]))),
        init_sigma2 = Sigma[j],
        init_delta2 = Delta[j],
        parameters_moves = parameters_moves,
        model_moves = model_moves,
        type_model_prior = type_model_prior,
        segmentation = segmentation, # if data segmentation should be applied
        segmentation_var = segmentation_var, # variable on which data is segmented
        n_segments = n_segments, # Fixed number of segments (if desired)
        segmentation_vector = segmentation_vector, # Fixed segmentation vector
        log = log,
        sigma2_in_betas_prior = sigma2_in_betas_prior
      )
      for (x in c('betas', 'sigma2s', 'delta2s', 'models', 'allocation_vectors')) {
        posterior_samples[[variables[j]]][[x]] <-  append(posterior_samples[[variables[j]]][[x]], res[[x]])
      }
    }
    posterior_samples_mat <- PosteriorSamplesMatrices(posterior_samples)
    M <- posterior_samples_mat$M
    B <- posterior_samples_mat$B
    Sigma <- posterior_samples_mat$Sigma
    Delta <- posterior_samples_mat$Delta
    missing_structure <- MissingFCDComponents(B, M, Sigma)
    #return(list(missing_mask, df, missing_structure,B, mean_0, sd_0))
    df <- MissingValuesUpdate(data = df, missing_mask = missing_mask, structure = missing_structure, mean_0 = mean_0, sd_0 = sd_0, sample_id = sample_id, time_id = time_id)
    if (epoch %% 100 == 0) {
      cat("Epoch:", (epoch * missing_update_period), "\n")
    }
    missing_values[[epoch]] <- as.matrix(df[,-c(1, 2)])[which(as.matrix(missing_mask[,-c(1, 2)]), arr.ind = TRUE)]
    #unlist(df[order(df[[sample_id]], df[[time_id]]), sort(setdiff(names(df), c(sample_id, time_id, "1")))][missing_mask[order(df[[sample_id]], df[[time_id]]), sort(setdiff(names(df), c(sample_id, time_id, "1")))] == TRUE], use.names = FALSE)
  }
  #posterior_samples
  list(posterior_samples, missing_values)
}

TemporalSmoothing <- function(data, sample_id = "Sample_Id", time_id = "Time", vars = NULL) {
  # Arrange data by sample and time
  data <- data %>% arrange(.data[[sample_id]], .data[[time_id]])
  
  if (is.null(vars)) {
    vars <- setdiff(names(data), c(sample_id, time_id))
  }
  
  imputed_data <- data
  
  # Step 1: Impute time 0 missing values using next available time point within sample
  time0_rows <- imputed_data[[time_id]] == 0
  
  for (var in vars) {
    # For each sample with missing at time 0
    samples_with_missing <- unique(imputed_data[[sample_id]][time0_rows & is.na(imputed_data[[var]])])
    
    for (s in samples_with_missing) {
      # Extract rows for this sample sorted by time
      sample_rows <- imputed_data[[sample_id]] == s
      sample_data <- imputed_data[sample_rows, c(time_id, var)]
      sample_data <- sample_data[order(sample_data[[time_id]]), ]
      
      # Find the earliest time > 0 with observed value for var
      next_val <- NA
      future_times <- sample_data[[time_id]] > 0
      
      if (any(future_times)) {
        future_vals <- sample_data[[var]][future_times]
        future_times_order <- sample_data[[time_id]][future_times]
        
        # Find first non-NA in future values
        idx <- which(!is.na(future_vals))[1]
        if (!is.na(idx)) {
          next_val <- future_vals[idx]
        }
      }
      
      # Impute missing at time 0 if next_val found
      if (!is.na(next_val)) {
        imputed_data[[var]][sample_rows & imputed_data[[time_id]] == 0 & is.na(imputed_data[[var]])] <- next_val
      } else {
        message(sprintf("Sample %s var %s: no future value to impute missing at time 0", s, var))
      }
    }
  }
  
  # Step 2: LOCF within each sample over time (fill forward)
  for (var in vars) {
    imputed_data <- imputed_data %>%
      group_by(.data[[sample_id]]) %>%
      arrange(.data[[time_id]]) %>%
      mutate(!!var := na.locf(.data[[var]], na.rm = FALSE)) %>%
      ungroup()
  }
  
  return(data.frame(imputed_data))
}

TemporalMice <- function(data, sample_id = "Sample_Id", time_id = "Time",
                                    vars = NULL, max_iter = 5) {
  data <- data[order(data[[sample_id]], data[[time_id]]), ]
  
  if (is.null(vars)) {
    vars <- setdiff(names(data), c(sample_id, time_id))
  }
  
  imputed_data <- data
  
  # Step 1: Initial mean imputation per (time, sample) group
  for (s in unique(imputed_data[[sample_id]])) {
    for (t in unique(imputed_data[[time_id]])) {
      idx <- imputed_data[[sample_id]] == s & imputed_data[[time_id]] == t
      for (var in vars) {
        vals <- imputed_data[[var]][idx]
        if (anyNA(vals)) {
          fallback_mean <- mean(imputed_data[[var]], na.rm = TRUE)
          imputed_data[[var]][idx & is.na(vals)] <- mean(vals, na.rm = TRUE)
          # If still NA (i.e. all missing at this (s,t)), use fallback
          imputed_data[[var]][idx & is.na(imputed_data[[var]])] <- fallback_mean
        }
      }
    }
  }
  
  # Step 2: Iterative modeling
  for (iter in seq_len(max_iter)) {
     #cat("Iteration", iter, "\n")
    
    for (target_var in vars) {
      #cat("  Imputing", target_var, "\n")
      
      for (s in unique(imputed_data[[sample_id]])) {
        sample_data <- imputed_data[imputed_data[[sample_id]] == s, ]
        time_vals <- sort(unique(sample_data[[time_id]]))
        
        for (i in 2:length(time_vals)) {
          t <- time_vals[i]
          t_prev <- time_vals[i - 1]
          
          # Data at time t and t-1
          data_t <- sample_data[sample_data[[time_id]] == t, ]
          data_t_prev <- sample_data[sample_data[[time_id]] == t_prev, ]
          
          # Skip if no data
          if (nrow(data_t) == 0 || nrow(data_t_prev) == 0) next
          
          # Create lagged covariate frame
          lagged_vars <- data_t_prev[, vars, drop = FALSE]
          colnames(lagged_vars) <- paste0(vars, "_lag")
          
          # Match by row index
          if (nrow(lagged_vars) != nrow(data_t)) {
            next  # skip if unequal size — can't match
          }
          
          df_model <- cbind(data_t[target_var], lagged_vars)
          names(df_model)[1] <- target_var
          
          complete_cases <- complete.cases(df_model)
          missing_mask <- is.na(data_t[[target_var]])
          
          if (sum(complete_cases) >= 2 && any(missing_mask)) {
            model <- lm(as.formula(paste(target_var, "~", paste0(colnames(lagged_vars), collapse = "+"))),
                        data = df_model[complete_cases, ])
            preds <- predict(model, newdata = lagged_vars[missing_mask, , drop = FALSE])
            
            # Update values
            idx_to_update <- which(imputed_data[[sample_id]] == s & imputed_data[[time_id]] == t)
            imputed_data[[target_var]][idx_to_update[missing_mask]] <- preds
          }
        }
      }
    }
  }
  
  return(imputed_data)
}

EvaluateMissingImputation <- function(data,
                                      df_with_missing,
                                      missing_posterior_samples,
                                      sample_id = "Sample_Id",
                                      time_id = "Time",
                                      burn_in = 500,
                                      credible_level = 0.95) {
  
  library(dplyr)
  library(ggplot2)
  
  # Ensure columns are sorted
  data <- data[, sort(names(data))]
  df_with_missing <- df_with_missing[, sort(names(df_with_missing))]
  
  # Create missing mask
  missing_mask <- df_with_missing %>%
    mutate(across(-c(sample_id, time_id), ~ is.na(.)))
  
  # Extract node (column) indices for missing values
  node_per_missing <- which(as.matrix(missing_mask[ , -c(1, 2)]), arr.ind = TRUE)[, 'col']
  
  # Extract true values at missing positions
  true_values <- as.matrix(data[ , -c(1, 2)])[which(as.matrix(missing_mask[ , -c(1, 2)]), arr.ind = TRUE)]
  
  # Extract posterior samples as a matrix: rows = iterations, cols = missing values
  posterior_matrix <- do.call(rbind, missing_posterior_samples)
  posterior_matrix <- posterior_matrix[burn_in:nrow(posterior_matrix), ]
  
  n_missings <- ncol(posterior_matrix)
  
  # Compute quantiles of true values in their posterior samples
  quantiles <- sapply(1:n_missings, function(i) {
    ecdf(posterior_matrix[ , i])(true_values[i])
  })
  
  # Compute credible intervals
  alpha <- (1 - credible_level) / 2
  credible_intervals <- t(sapply(1:n_missings, function(i) {
    quantile(posterior_matrix[ , i], probs = c(alpha, 1 - alpha))
  }))
  
  # Count how many missing per variable to assign a local index
  var_names <- colnames(data)[-c(1, 2)]
  variable_names <- factor(var_names[node_per_missing], levels = var_names)
  variable_counts <- table(variable_names)
  
  # Create local index within each variable
  local_indices <- unlist(lapply(variable_counts, function(n) seq_len(n)))
  
  # Combine into data frame for plotting
  distributions_df <- data.frame(
    Variable = variable_names,
    LocalIndex = local_indices,
    Lower = credible_intervals[ , 1],
    Upper = credible_intervals[ , 2],
    true_values = true_values
  )
  
  # Plot with variable on x-axis and some dodging between points
  p <- ggplot(distributions_df, aes(x = Variable, group = interaction(Variable, LocalIndex))) +
    geom_linerange(aes(ymin = Lower, ymax = Upper), 
                   position = position_dodge(width = 0.6), color = "gray70") +
    geom_point(aes(y = true_values, color = Variable), 
               position = position_dodge(width = 0.6), size = 1.5) +
    labs(title = paste0("Posterior ", round(100 * credible_level), "% Credible Intervals vs. True Values"),
         x = "Variable",
         y = "Value",
         color = "Variable") +
    theme_minimal(base_size = 14)
  
  # Coverage statistic
  coverage <- mean(quantiles >= alpha & quantiles <= (1 - alpha))
  
  list(
    plot = p,
    coverage = coverage,
    quantiles = quantiles,
    distributions_df = distributions_df
  )
}


# Function to compute PSRF for each missing value
ComputePSRFCurve <- function(missing_matrices, threshold = 1.1, start_iter = 10, step = 1) {
  n_iter <- nrow(missing_matrices[[1]])
  n_missings <- ncol(missing_matrices[[1]])
  iter_points <- seq(start_iter, n_iter, by = step)
  
  psrf_over_time <- numeric(length(iter_points))
  
  for (t_idx in seq_along(iter_points)) {
    t <- iter_points[t_idx]
    
    # Compute PSRF at iteration t for each missing value
    psrf_values <- sapply(1:n_missings, function(j) {
      mcmc_list <- mcmc.list(lapply(missing_matrices, function(mat) mcmc(mat[1:t, j])))
      gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1]
    })
    
    # Fraction of PSRFs < threshold
    psrf_over_time[t_idx] <- mean(psrf_values < threshold)
  }
  
  # Return data frame for plotting
  data.frame(
    Iteration = iter_points,
    FractionBelowThreshold = psrf_over_time
  )
}
