StandardizeData <- function(df, cols = colnames(df)[!colnames(df) %in% c('1','Sample_Id', 'Time')]){
  df[,cols] <- data.frame(scale(df[,cols]))
  df
}

StandardizeBeta <- function(original_betas, original_data) {
  # Get the number of variables (predictors) including the intercept
  n_vars <- ncol(original_data) - 1  # Exclude the intercept column ('1')
  
  # Calculate the mean and standard deviation for each predictor (excluding 'Y')
  mu_original <- colMeans(original_data[, -ncol(original_data)])  # Exclude 'Y' column
  sigma_original <- apply(original_data[, -ncol(original_data)], 2, sd)  # Exclude 'Y' column
  # Standard deviation of the outcome variable Y
  mu_Y <- mean(original_data$Y)
  sigma_Y <- sd(original_data$Y)
  
  # Transform the betas for the predictors (excluding intercept)
  standardized_betas <- original_betas * sigma_original / sigma_Y  # Adjusting betas for X_j's
  
  # Adjust intercept: the intercept is affected by the means of the predictors, their SDs, and the SD of Y
  
  standardized_betas[1] <- (original_betas[1] - mu_Y + sum(original_betas[2:n_vars] * mu_original[2:n_vars])) / sigma_Y 
  standardized_betas
}

MC_selection = function(betas, sigma2s, delta2s, burn_in = floor(length(betas) / 2), thin_out = 5){
  idx = seq(burn_in + 1, length(betas), by = thin_out)  # start from burn_in + 1
  list(
    betas = betas[idx],
    sigma2s = sigma2s[idx],
    delta2s = delta2s[idx]
  )
}

RJ_MC_selection = function(betas, sigma2s, delta2s, models, allocation_vectors, covariates = NULL, burn_in = round(length(betas) / 2), thin_out = 5){
  #Return thinned samples post-burn-in phase
  if (length(covariates) == 0){
    list(betas = betas[seq(burn_in, length(betas), by = thin_out)], sigma2s =
           sigma2s[seq(burn_in, length(sigma2s), by = thin_out)], delta2s = delta2s[seq(burn_in, length(delta2s), by = thin_out)], models = models[seq(burn_in, length(models), by = thin_out)], allocation_vectors = allocation_vectors[seq(burn_in, length(models), by = thin_out)])
  }
  else{
    list(betas = betas[seq(burn_in, length(betas), by = thin_out)], sigma2s =
           sigma2s[seq(burn_in, length(sigma2s), by = thin_out)], delta2s = delta2s[seq(burn_in, length(delta2s), by = thin_out)], models = models[seq(burn_in, length(models), by = thin_out)], allocation_vectors = allocation_vectors[seq(burn_in, length(models), by = thin_out)], covariates = covariates)
    
  }
}

CovariateConfigurations <-
  function(covariates, include_empty = FALSE) {
    k <- length(covariates)
    
    # Generate each combination of the covariates in the model
    combos <- lapply(0:k, function(i) {
      combn(covariates, i, simplify = FALSE)
    })
    all_configs <- unlist(combos, recursive = FALSE)
    
    # If include_empty also the model with just the Intercept is considered
    if (!include_empty) {
      all_configs <- Filter(length, all_configs)
    }
    all_configs
  }

DIC <-
  function(mcmc_output,
           data,
           outcome = "Y",
           covariates = colnames(data)[!colnames(data) %in% outcome]) {
    # Extract posteriors
    beta_samples <- mcmc_output$betas
    sigma2_samples <- mcmc_output$sigma2s
    
    # Prepare design matrix
    X <- as.matrix(data[, covariates])
    Y <- data[[outcome]]
    N <- length(Y)
    S <- length(beta_samples)  # number of posterior samples
    
    # Compute deviances for each posterior sample
    deviances <- numeric(S)
    for (s in 1:S) {
      beta <- beta_samples[[s]]
      sigma2 <- sigma2_samples[[s]]
      resid <- Y - X %*% beta
      deviances[s] <- N * log(2 * pi * sigma2) + sum(resid ^ 2) / sigma2
    }
    
    # Posterior means
    beta_mean <- Reduce("+", beta_samples) / S
    sigma2_mean <- mean(unlist(sigma2_samples))
    resid_mean <- Y - X %*% beta_mean
    dev_at_mean <-
      N * log(2 * pi * sigma2_mean) + sum(resid_mean ^ 2) / sigma2_mean
    
    # DIC Computation
    DIC <- 2 * mean(deviances) - dev_at_mean
    cat(DIC)
    cat('\n\n')
    DIC
}

InclusionProbsToTransitionNetwork <- function(inclusion_probs_list, threshold = 0.05, default_variance = 1) {
  lapply(names(inclusion_probs_list), function(target_var) {
    probs <- inclusion_probs_list[[target_var]]
    
    # Filter based on threshold
    selected <- probs[abs(probs) > threshold]
    
    # Clean parent names (remove '_t-1'), keep '1' untouched
    cleaned_parents <- gsub("_t-1$", "", names(selected))
    
    # Clean target variable name (remove '_t')
    cleaned_target <- gsub("_t$", "", target_var)
    
    list(
      var_name = cleaned_target,
      parents = cleaned_parents,
      params = as.numeric(selected),
      variance = default_variance
    )
  }) |>
    setNames(gsub("_t$", "", names(inclusion_probs_list)))
}

Summarize_DBNSamples <- function(DBN_results, 
                                 true_values_list = NULL,
                                 varnames = NULL,
                                 burn_in_rate = 0.5,
                                 thin_out = 5,
                                 ncols = 3) {
  
  summaries <- list()
  
  for (j in seq_along(DBN_results)) {
    result <- DBN_results[[j]]
    varname <- if (!is.null(varnames)) varnames[j] else paste0("X_", j)
    
    # Apply MCMC selection (burn-in and thinning)
    selected <- RJ_MC_selection(
      betas = result$betas,
      sigma2s = result$sigma2s,
      delta2s = result$delta2s,
      models = result$models,
      allocation_vectors = result$allocation_vectors,
      burn_in = round(length(result$betas) * burn_in_rate),
      thin_out = thin_out
    )
    # Remove empty model if present in the thinned sample
    beta_samples <- selected$betas
    model_samples <- selected$betas
    # Posterior plot for betas
    cat("Variable:", varname, "\n")
    beta_summary <- PosteriorDistribution(
      posterior_samples = beta_samples,
      true_values = if (!is.null(true_values_list)) true_values_list[[j]] else NULL,
      ncols = ncols,
      covariates = result$covariates
    )
    
    # Posterior plot for sigma² and delta²
    sigma_delta_summary <- PosteriorDistributionSigmadelta(
      mcmc_output = selected,
      true_sigma2 = NULL,
      true_delta2 = NULL
    )
    
    summaries[[varname]] <- list(
      beta_summary = beta_summary,
      sigma_delta_summary = sigma_delta_summary
    )
  }
  
  return(summaries)
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