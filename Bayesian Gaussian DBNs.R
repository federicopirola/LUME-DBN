library(dplyr)
library(ggplot2)
library(ggforce)
library(grid)
library(latex2exp)
library(ggtext)

TransitionNetwork <- function(k,
                              max_parent_set = 3,
                              a1 = 0.1, 
                              a2 = 1,
                              V = 1,
                              var_set = sapply(1:k, function(i) {paste0("X_", i)})){
  topological_order <- sample(var_set)
  transition_network <- list()
  for (i in 1:k){
    if (i == 1){
      transition_network[[topological_order[i]]] <- list(var_name = topological_order[i],
                                                         parents = c('1'),
                                                         params = sample(c(-1, 1), 1, replace = TRUE) * runif(1, a1, a2),
                                                         variance = V)
    }
    else{
      n_parents <- sample(0:min(i-1,max_parent_set),1)
      parent_set <- topological_order[sample(1:(i-1), n_parents)]
      transition_network[[topological_order[i]]] <- list(var_name = topological_order[i],
                                                         parents = c('1',parent_set),
                                                         params = sample(c(-1, 1), n_parents + 1, replace = TRUE) * runif(n_parents + 1, a1, a2),
                                                         variance = V)
    }
  }
  transition_network[order(names(transition_network))]
}

Network0 <- function(k,
                     a1 = 0.1,
                     a2 = 1,
                     V = 1,
                     var_set = sapply(1:k, function(i) {paste0("X_", i)})){
  setNames(
    lapply(var_set, function(i) list(
      variable = i,
      params = sample(c(-1, 1), 1) * runif(1, a1, a2),
      variance = V
    )),
    var_set
  )
}

DBNGeneration <- function(k,
                          max_parent_set = 3,
                          a1 = 0.1, 
                          a2 = 1,
                          V = 1,
                          var_set = sapply(1:k, function(i) {paste0("X_", i)})){
  list(Network0 = Network0(k = k, a1 = a1, a2 = a2, V = V, var_set = var_set),
       TransitionNetwork = TransitionNetwork(k = k, max_parent_set = max_parent_set, a1 = a1, a2 = a2, V = V, var_set = var_set))
}

DataGenerationFromDBN <- function(N, T, dbn, scaling = FALSE) {
  initial_net <- dbn$Network0
  transition_net <- dbn$TransitionNetwork
  
  k <- length(initial_net)
  variable_names <- names(initial_net)
  
  # --- Time 0: Generate from Network0 ---
  initial_data <- replicate(N, {
    sapply(initial_net, function(model) {
      rnorm(1, mean = model$params, sd = sqrt(model$variance))
    })
  })
  initial_df <- data.frame(Sample_Id = 1:N, Time = 0, t(initial_data))
  colnames(initial_df)[3:ncol(initial_df)] <- variable_names
  
  df <- initial_df
  
  # --- Time 1 to T: Generate from TransitionNetwork ---
  for (t in 1:T) {
    past_df <- df[df$Time == (t - 1), ]
    new_df <- data.frame(Sample_Id = 1:N, Time = t)
    
    for (var in variable_names) {
      model <- transition_net[[var]]
      parents <- model$parents
      params <- model$params
      
      intercept <- params[1]
      parent_vars <- parents[-1]
      parent_params <- if (length(parents) > 1) params[-1] else numeric(0)
      
      if (length(parent_vars) > 0) {
        predictors <- past_df[, parent_vars, drop = FALSE]
        if (scaling) {
          predictors <- scale(predictors)
        }
        lin_comb <- intercept + as.numeric(as.matrix(predictors) %*% parent_params)
      } else {
        lin_comb <- rep(intercept, N)
      }
      
      new_df[[var]] <- rnorm(N, mean = lin_comb, sd = sqrt(model$variance))
    }
    
    df <- rbind(df, new_df)
  }
  
  return(df)
}

NH_DBNetwork <- function(k,
                         T,
                         changepoints = NULL,
                         max_parent_set = 3,
                         a1 = 0.1, 
                         a2 = 1,
                         V = 1,
                         var_set = sapply(1:k, function(i) paste0("X_", i)),
                         param_type = c("independent", "same_sign", "flip_sign")) {
  
  param_type <- match.arg(param_type)
  
  # -- Handle changepoints: time t in changepoints means new params start at t
  if (is.null(changepoints)) {
    n_segments <- 3
    changepoints <- round(seq(1, T, length.out = n_segments + 1))[-1]
  }
  changepoints <- sort(unique(changepoints))
  segment_starts <- c(1, changepoints)
  segment_ends <- c(changepoints - 1, T)
  n_segments <- length(segment_starts)
  
  # Fixed structure
  topological_order <- sample(var_set)
  structure <- list()
  for (i in 1:k) {
    if (i == 1) {
      structure[[topological_order[i]]] <- c("1")
    } else {
      n_parents <- sample(0:min(i - 1, max_parent_set), 1)
      structure[[topological_order[i]]] <- c("1", topological_order[sample(1:(i - 1), n_parents)])
    }
  }
  
  # Parameter sets
  parameter_sets <- vector("list", n_segments)
  for (s in 1:n_segments) {
    parameter_sets[[s]] <- list()
    for (var in topological_order) {
      parents <- structure[[var]]
      n_parents <- length(parents) - 1
      
      if (s == 1 || param_type == "independent") {
        params <- sample(c(-1, 1), n_parents + 1, replace = TRUE) * runif(n_parents + 1, a1, a2)
      } else {
        prev_params <- parameter_sets[[s - 1]][[var]]$params
        if (param_type == "same_sign") {
          signs <- sign(prev_params)
          params <- signs * runif(n_parents + 1, a1, a2)
        } else if (param_type == "flip_sign") {
          params <- -prev_params
        }
      }
      
      parameter_sets[[s]][[var]] <- list(
        var_name = var,
        parents = parents,
        params = params,
        variance = V
      )
    }
  }
  
  list(
    Network0 = Network0(k = k, a1 = a1, a2 = a2, V = V, var_set = var_set),
    changepoints = changepoints,
    segment_starts = segment_starts,
    TransitionNetworks = parameter_sets
  )
}


DataGenerationFromNHDBN <- function(N, T, nh_dbn, scaling = FALSE) {
  initial_net <- nh_dbn$Network0
  variable_names <- names(initial_net)
  segment_starts <- nh_dbn$segment_starts
  segments <- nh_dbn$TransitionNetworks
  n_segments <- length(segments)
  
  # Time 0
  initial_data <- replicate(N, {
    sapply(initial_net, function(model) {
      rnorm(1, mean = model$params, sd = sqrt(model$variance))
    })
  })
  df <- data.frame(Sample_Id = 1:N, Time = 0, t(initial_data))
  colnames(df)[3:ncol(df)] <- variable_names
  
  # Time 1 to T
  for (t in 1:T) {
    past_df <- df[df$Time == (t - 1), ]
    new_df <- data.frame(Sample_Id = 1:N, Time = t)
    
    # Identify segment
    s_idx <- max(which(segment_starts <= t))
    current_params <- segments[[s_idx]]
    
    for (var in variable_names) {
      model <- current_params[[var]]
      parents <- model$parents
      params <- model$params
      intercept <- params[1]
      parent_vars <- parents[-1]
      parent_params <- if (length(parents) > 1) params[-1] else numeric(0)
      
      if (length(parent_vars) > 0) {
        predictors <- past_df[, parent_vars, drop = FALSE]
        if (scaling) {
          predictors <- scale(predictors)
        }
        lin_comb <- intercept + as.numeric(as.matrix(predictors) %*% parent_params)
      } else {
        lin_comb <- rep(intercept, N)
      }
      
      new_df[[var]] <- rnorm(N, mean = lin_comb, sd = sqrt(model$variance))
    }
    
    df <- rbind(df, new_df)
  }
  
  return(df)
}




LaggedDatasetGeneration <- function(df, sample_id = 'Sample_Id', time_id = 'Time',lag = 1){
  # Get max time from data
  T_max <- max(df$Time)
  # Get variable names (assumes all columns except Sample_Id and Time are variables)
  variables <- setdiff(names(df), c(sample_id, time_id))
  k <- length(variables)
  
  # Create lagged dataset for DBN
  lagged_data <- do.call(rbind, lapply(lag:T_max, function(t) {
    current <- df[df$Time == t, c(sample_id, time_id, variables)]
    colnames(current) <- c(sample_id, time_id,paste0(variables, "_t"))
    lagged_df <- current
    for (time_lag in 1:lag) {
      temp <- df[df$Time == (t - time_lag), variables]
      colnames(temp) <- paste0(variables, "_t-", time_lag)
      lagged_df <- cbind(lagged_df, temp)
    }
    lagged_df
  }))  # prepend intercept column named "1" 
 lagged_data
}

LaggedDatasetGeneration2 <- function(df, sample_id = "Sample Id", time_id = "Time", lag = 1) {
  variables <- setdiff(names(df), c(sample_id, time_id))
  
  # Result storage
  result_list <- list()
  
  # Unique samples
  sample_vals <- unique(df[[sample_id]])
  
  for (s in sample_vals) {
    this_sample <- df[df[[sample_id]] == s, ]
    times <- this_sample[[time_id]]
    
    valid_times <- times[times %in% (lag + min(times)):max(times)]
    
    for (t in valid_times) {
      current_row <- this_sample[this_sample[[time_id]] == t, ]
      # Skip if current row or variables have NA or row missing
      if (nrow(current_row) == 0 || any(is.na(current_row[, variables, drop = FALSE]))) next
      
      lagged_df <- current_row[, c(sample_id, time_id), drop = FALSE]
      
      # Current variables with _t suffix
      current_vars <- current_row[, variables, drop = FALSE]
      colnames(current_vars) <- paste0(variables, "_t")
      lagged_df <- cbind(lagged_df, current_vars)
      
      # Add lagged variables
      missing_lag <- FALSE
      for (l in 1:lag) {
        lagged_row <- this_sample[this_sample[[time_id]] == (t - l), variables, drop = FALSE]
        if (nrow(lagged_row) == 0 || any(is.na(lagged_row))) {
          missing_lag <- TRUE
          break
        }
        colnames(lagged_row) <- paste0(variables, "_t-", l)
        lagged_df <- cbind(lagged_df, lagged_row)
      }
      
      if (!missing_lag) {
        result_list[[length(result_list) + 1]] <- lagged_df
      }
    }
  }
  
  # Combine all rows
  if (length(result_list) > 0) {
    result <- do.call(rbind, result_list)
    rownames(result) <- NULL
    return(result)
  } else {
    return(NULL)
  }
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



# TO DO: Init Model Matrix - Init Betas Matrix - Init Allocation List - Init Sigma2 Vector - Betas Mean matrix
           
RJMCMC_DBN <- function(df,
                       sample_id = 'Sample_Id',
                       time_id = 'Time',
                       parameters_moves = 'Gibbs',
                       model_moves = 'MH-Marginal',
                       type_model_prior = 'Uniform',
                       type_segment_prior = 'Uniform_card', # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
                       type_coupling = 'piece-wise', # type of parameter coupling ('piece-wise' or 'global')
                       coupling = FALSE, # if parameter coupling should be applied
                       standardized = FALSE,
                       lagged = TRUE,
                       time_lags = 1,
                       self_loop = FALSE,
                       p_MCMC = 1, # probability of doing an MCMC move
                       p_RJMCMC = 0.1,
                       p_ChangePoint = 0.1, # probability of doing a Change Point Detection move
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
                       init_betas = rep(0, (length(covars))),
                       init_xi2 = 1, # initial xi2
                       init_allocation_vector = NaN, # initial segment set
                       segmentation = FALSE, # if data segmentation should be applied
                       segmentation_var = 'Time', # variable on which data is segmented
                       n_segments = NaN, # Fixed number of segments (if desired)
                       segmentation_vector = NaN, # Fixed segmentation vector
                       log = TRUE,
                       sigma2_in_betas_prior = FALSE
                       ) {
  # Get variable names (assumes all columns except Sample_Id and Time are variables)
  if (lagged){
    variables <- setdiff(names(df), c(sample_id, time_id))
  }
  else{
    variables <- gsub('_t','',names(df)[grepl("_t$", names(df))])
  }
  k <- length(variables)
  if (standardized){
    df[,variables] <- StandardizeData(df[,variables])
  }
  # Create lagged dataset for DBN
  if (lagged){
    lagged_data <- LaggedDatasetGeneration(df, lag = time_lags, sample_id = sample_id, time_id = time_id)
    lagged_data <- cbind(lagged_data[, 1:2], '1' = 1, lagged_data[, 3:ncol(lagged_data)])
  }
  else{
    lagged_data <- cbind(df[, 1:2], '1' = 1, df[, 3:ncol(df)])
  }
  # Fit one model for each variable at time t
  posterior_samples <- list()
  for (j in 1:k) {
    outcome_var <- paste0(variables[j], "_t")
    if (!self_loop) {
      covars <- unlist(lapply(1:time_lags, function(lag) paste0(variables[-j], "_t-", lag)))
    }
    else{
      covars <- unlist(lapply(1:time_lags, function(lag) paste0(variables, "_t-", lag)))  
    }
    # Define covariates as all variables at t-1
    cat(outcome_var,'\n\n')
    # Include intercept if required
    covars <- c("1", covars)
    res <- RJMCMC(
      data = lagged_data[,c(sample_id, time_id, covars, outcome_var)],
      epochs = epochs,
      v = v,
      type_segment_prior = type_segment_prior, # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
      type_coupling = type_coupling, # type of parameter coupling ('piece-wise' or 'global')
      coupling = coupling, # if parameter coupling should be applied
      p_MCMC = p_MCMC, # probability of doing an MCMC move
      p_RJMCMC = p_RJMCMC,
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
      covariates = covars,
      init_model = character(0),
      init_betas = init_betas,
      betas_mean = rep(0, (length(covars))),
      init_sigma2 = init_sigma2,
      init_delta2 = init_delta2,
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
    posterior_samples[[gsub("_t", "", outcome_var)]] <- res
    posterior_samples[[gsub("_t", "", outcome_var)]][['covariates']] <- gsub("_t-1", "", covars)
  }
  
  posterior_samples
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

TransitionNetworkPlot <- function(transition_net,
                                  threshold = 0.05,
                                  node_size = 10,
                                  curvature = 0,
                                  text_size = 4.5) {
  library(ggplot2)
  library(geomtextpath)
  
  vars <- names(transition_net)
  k <- length(vars)
  
  # Layouts for t-1 and t
  layout_tminus1 <- data.frame(
    var = vars,
    label = vars,
    time = "t-1",
    x = 0,
    y = seq(k, 1)
  )
  layout_t <- data.frame(
    var = vars,
    label = vars,
    time = "t",
    x = 3,  # wider spacing between time steps
    y = seq(k, 1)
  )
  nodes <- rbind(layout_tminus1, layout_t)
  
  # Create edges (no labels)
  edges <- do.call(rbind, lapply(transition_net, function(entry) {
    target <- entry$var_name
    parents <- entry$parents
    params <- entry$params
    data.frame(
      from = parents,
      to = target,
      x = layout_tminus1$x[match(parents, layout_tminus1$var)],
      y = layout_tminus1$y[match(parents, layout_tminus1$var)],
      xend = layout_t$x[match(target, layout_t$var)],
      yend = layout_t$y[match(target, layout_t$var)],
      strength = params
    )
  }))
  
  # Filter by threshold
  edges <- subset(edges, abs(strength) > 0)
  edges$lty <- ifelse(abs(edges$strength) >= threshold, "solid", "dotted")
  
  # Plot
  ggplot() +
    # Edges
    geom_curve(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend, linetype = lty),
      curvature = curvature,
      arrow = arrow(length = unit(0.2, "cm")),
      color = "gray30",
      show.legend = FALSE
    ) +
    # Time slice borders
    annotate("rect", xmin = -0.6, xmax = 0.6, ymin = 0.5, ymax = k + 0.5,
             fill = NA, color = "black", size = 0.6) +
    annotate("rect", xmin = 2.4, xmax = 3.6, ymin = 0.5, ymax = k + 0.5,
             fill = NA, color = "black", size = 0.6) +
    # Time labels above boxes
    annotate("text", x = 0, y = k + 1.2, label = expression(t - 1), size = 6) +
    annotate("text", x = 3, y = k + 1.2, label = expression(t), size = 6) +
    # Nodes
    geom_point(data = nodes, aes(x = x, y = y, fill = time),
               size = node_size, shape = 21, color = "black") +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              size = text_size) +
    scale_fill_manual(values = c("t-1" = "skyblue", "t" = "tomato")) +
    coord_fixed(xlim = c(-1, 4), ylim = c(0.5, k + 1.5)) +
    theme_void() +
    ggtitle("Transition Network") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}