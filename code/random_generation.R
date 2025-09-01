CompleteLinearRegressionModelGeneration = function(n_strong_predictors = 3,
                                                   n_week_predictors = 0,
                                                   n_non_predictors = 0,
                                                   a_1 = 3,
                                                   a_2 = 5,
                                                   b_1 = 0.2,
                                                   b_2 = 0.5,
                                                   M = rep(0, n_strong_predictors + n_week_predictors + n_non_predictors),
                                                   V = rep(1, n_strong_predictors + n_week_predictors + n_non_predictors)) {
  # Randomly generating Betas for strong predictors from the Uniform distribution in the space [-a_2,-a_1] U [a_1,a_2]
  betas <-
    ifelse(rbinom(n_strong_predictors, 1, 0.5) == 1, 1,-1) * runif(n_strong_predictors, min = a_1, max = a_2)
  
  # Randomly generating Betas for week predictors from the Uniform distribution in the space [-b_2,-b_1] U [b_1,b_2]
  betas <- c(betas, ifelse(rbinom(n_week_predictors, 1, 0.5) == 1, 1,-1) * runif(n_week_predictors, min = b_1, max = b_2))
  
  # Setting betas for non predictors equal to 0  
  betas <- c(betas, rep(0, n_non_predictors))
  
  # Shuffle the parameters
  betas <- sample(betas)
  
  # Saving all the parameters of the model in a list
  list(
    # Intercept + Linear Coefficients
    betas = betas,
    # Means
    means = M,
    # Standard Deviations
    variances = V
  )
}

DataGeneration = function(N = 100, parameters) {
  # Number of Variables
  p <- length(parameters$means)
  if (length(parameters$variances) != p) {
    stop("Means and variances must have the same length.")
  }
  
  # Generating Covariates Values given their means and variances
  data <-
    mapply(function(mean, variance)
      rnorm(N, mean = mean, sd = sqrt(variance)),
      parameters$means[1:(p - 1)],
      parameters$variances[1:(p - 1)])
  
  df <- as.data.frame(cbind(1,data))
  colnames(df) <- c(1,paste0("X", 1:(p - 1)))
  
  # Generating Outcome Values by combining normals: linear combination of the covariates plus a normal error
  Y <-
    as.matrix(df) %*% parameters$betas + rnorm(N, parameters$means[p], sqrt(parameters$variances[p]))
  df$Y <- Y
  return(df)
}

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
