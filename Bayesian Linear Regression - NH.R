library(emdbook)
library(extraDistr)
library(ggplot2)
library(geomtextpath)
library(reshape2)
library(MASS)

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

ModelPlot <- function(strengths,
                      covariates = c('1',paste0("X", 1:(length(strengths)-1))),
                      outcome = 'Y',
                      threshold = 0.2) {
  k = length(strengths)
  # Defining covariates names
  covariates_df <- data.frame(name = covariates,
                              x = 0,
                              y = seq(k, 1))
  
  # Defining outcome node
  y_node <- data.frame(name = outcome,
                       x = 2 * k / 3,
                       y = mean(covariates_df$y))
  
  nodes <- rbind(covariates_df, y_node)
  
  # Create edges only where strength != 0
  edges <- data.frame(
    x = covariates_df$x,
    y = covariates_df$y,
    xend = y_node$x,
    yend = y_node$y,
    label = round(strengths, 2),
    strength = strengths
  )
  
  # Filter out zero-strength edges
  edges <- subset(edges, strength != 0)
  
  # Determine line type
  edges$lty <- ifelse(abs(edges$strength) >= threshold, "dotted", "solid")
  
  # Plot
  ggplot() +
    # Curved arrows with labels
    geom_textcurve(
      data = edges,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        label = label,
        linetype = lty
      ),
      curvature = 0,
      arrow = arrow(length = unit(0.2, "cm")),
      size = 4,
      color = "gray20",
      text_smoothing = 30,
      show.legend = FALSE
    ) +
    # Nodes
    geom_point(
      data = nodes,
      aes(x = x, y = y, color = name),
      size = 16,
      show.legend = FALSE
    ) +
    # Node labels
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = name),
      color = "white",
      size = 5.5
    ) +
    scale_color_manual(values = c(rep("darkgreen", k), "tomato")) +
    theme_void() +
    coord_fixed(xlim = c(-0.5, 0.5 + 2 * k / 3),
                ylim = c(0.5, k + 0.5))
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

PosteriorDistribution = function(posterior_samples,
                                 true_values,
                                 ncols = 3,
                                 covariates = NULL) {
  covariates[1] = 'Intercept'
  
  nparams = length(posterior_samples[[1]])
  nrows = ceiling(nparams / ncols)
  
  # Combine all MCMC chains into one matrix
  mcmc_list <-
    lapply(posterior_samples, function(x)
      matrix(x, ncol = nparams, byrow = TRUE))
  mcmc_matrix <- do.call(rbind, mcmc_list)
  
  # Compute posterior means and 95% credible intervals
  posterior_means <- colMeans(mcmc_matrix)
  credible_intervals <-
    apply(mcmc_matrix, 2, function(x)
      quantile(x, probs = c(0.025, 0.975)))
  colnames(credible_intervals) <- covariates
  
  # Set plotting layout
  par(mfrow = c(nrows, ncols))
  print(credible_intervals)
  covariates
  
  # Iteration over Beta_j
  for (j in 1:nparams) {
    samples <- mcmc_matrix[, j]
    
    if (!is.null(true_values)) {
      # Definying grid space based on random values and true values
      x_min <- min(min(samples), true_values[[j]]) - 0.05
      x_max <- max(max(samples), true_values[[j]]) + 0.05
      
      # Plotting Histogram
      hist_info <- hist(
        samples,
        main = bquote("Posterior of " ~ beta[.(j - 1)]),
        xlab = bquote(beta[.(j - 1)]),
        col = "gray",
        border = "white",
        breaks = 20,
        freq = TRUE,
        xlim = c(x_min, x_max)
      )
      
      # Adding a Vertical Line for the true observed value
      abline(v = true_values[[j]],
             col = "green",
             lwd = 2)
      y_max <- max(hist_info$counts)
    }
    else{
      hist_info <- hist(
        samples,
        main = bquote("Posterior of " ~ beta[.(j - 1)]  (.(covariates[j]))),
        xlab = bquote(beta[.(j - 1)]),
        col = "gray",
        border = "white",
        breaks = 20,
        freq = TRUE
      )
      y_max <- max(hist_info$counts)
    }
    
    
    # Adding a Vertical Line for the Posterior mean
    abline(
      v = posterior_means[[j]],
      col = "red",
      lwd = 2,
      lty = 2
    )
    
    # Adding a Shaded Rectangle for the Credible Interval based on MCMC simulations
    rect(
      xleft = credible_intervals[1, j],
      xright = credible_intervals[2, j],
      ybottom = 0,
      ytop = y_max,
      col = rgb(0, 0, 0.2, alpha = 0.2),
      border = NA
    )
  }
  
  par(mfrow = c(1, 1))
  credible_intervals
}

PosteriorDistributionSegmented <- function(posterior_samples,
                                           true_values = NULL,
                                           ncols = 3,
                                           covariates = NULL,
                                           segment_labels = NULL,
                                           colors = NULL,
                                           show_mean = FALSE,
                                           n_changepoints = NULL) {
  # Group posterior_samples by number of changepoints (segments - 1)
  segment_counts <- sapply(posterior_samples, length)
  changepoint_counts <- segment_counts - 1
  
  if (!is.null(n_changepoints)) {
    # Filter by fixed number of changepoints
    keep <- which(changepoint_counts == n_changepoints)
    if (length(keep) == 0) stop("No samples with the specified number of changepoints.")
    posterior_samples <- posterior_samples[keep]
    changepoint_counts <- changepoint_counts[keep]
    n_changepoints_list <- n_changepoints
  } else {
    n_changepoints_list <- sort(unique(changepoint_counts))
  }
  
  all_summaries <- list()
  
  for (k in n_changepoints_list) {
    samples_k <- posterior_samples[changepoint_counts == k]
    summary_k <- .plot_segmented_distributions(samples_k, true_values, ncols, covariates,
                                               segment_labels = paste("Segment", 1:(k+1)),
                                               colors = rainbow(k+1, alpha = 0.7),
                                               show_mean = show_mean)
    all_summaries[[paste0(k, "_changepoints")]] <- summary_k
  }
  
  return(all_summaries)
}

.plot_segmented_distributions <- function(posterior_samples, true_values, ncols, covariates,
                                          segment_labels, colors, show_mean) {
  n_segments <- length(posterior_samples[[1]])
  n_iters <- length(posterior_samples)
  n_params <- length(posterior_samples[[1]][[1]])
  
  if (is.null(covariates)) covariates <- paste0("V", 1:n_params)
  covariates[1] <- "Intercept"
  
  # Convert to matrix per segment
  segment_samples <- vector("list", n_segments)
  for (s in 1:n_segments) {
    segment_samples[[s]] <- matrix(NA, nrow = n_iters, ncol = n_params)
  }
  for (iter in seq_along(posterior_samples)) {
    for (s in 1:n_segments) {
      segment_samples[[s]][iter, ] <- as.numeric(posterior_samples[[iter]][[s]])
    }
  }
  
  nrows <- ceiling(n_params / ncols)
  par(mfrow = c(nrows, ncols))
  
  summary_list <- list()
  
  for (j in 1:n_params) {
    densities <- lapply(segment_samples, function(mat) density(mat[, j]))
    x_min <- min(sapply(densities, function(d) min(d$x)))
    x_max <- max(sapply(densities, function(d) max(d$x)))
    y_max <- max(sapply(densities, function(d) max(d$y)))
    
    # Set plot limits including space for boxplots
    box_gap <- 0.06 * y_max
    y_box_base <- - (length(segment_samples) + 1) * box_gap
    y_lim <- c(y_box_base, 1.1 * y_max)
    
    # Initialize plot
    plot(0, 0,
         type = "n",
         xlim = c(x_min, x_max),
         ylim = y_lim,
         main = bquote("Posterior of " ~ beta[.(j - 1)] ~ "(" ~ .(covariates[j]) ~ ")"),
         xlab = bquote(beta[.(j - 1)]),
         ylab = "Density")
    
    # Plot densities
    for (s in 1:n_segments) {
      samples_j <- segment_samples[[s]][, j]
      dens <- density(samples_j)
      lines(dens, col = colors[s], lwd = 2)
    }
    
    # Add boxplots stacked vertically
    box_height <- 0.025 * y_max
    
    param_summary <- data.frame(
      Segment = segment_labels,
      Mean = NA,
      Median = NA,
      CI_2.5 = NA,
      CI_97.5 = NA
    )
    
    for (s in 1:n_segments) {
      samples_j <- segment_samples[[s]][, j]
      q <- quantile(samples_j, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
      stat_value <- if (show_mean) mean(samples_j) else q[3]
      y_box <- y_box_base + s * box_gap
      
      # Whiskers
      segments(q[1], y_box, q[5], y_box, col = colors[s], lwd = 2)
      # Box (IQR)
      rect(q[2], y_box - box_height / 2,
           q[4], y_box + box_height / 2,
           col = adjustcolor(colors[s], alpha.f = 0.5),
           border = colors[s])
      # Median/mean point
      points(stat_value, y_box, pch = 19, col = colors[s], cex = 1.2)
      # Segment label
      text(x_min, y_box, labels = segment_labels[s], pos = 2, cex = 0.8)
      
      # Save summaries
      param_summary$Mean[s] <- mean(samples_j)
      param_summary$Median[s] <- q[3]
      param_summary$CI_2.5[s] <- q[1]
      param_summary$CI_97.5[s] <- q[5]
    }
    
    if (!is.null(true_values)) {
      abline(v = true_values[[j]], col = "green", lwd = 2)
    }
    
    summary_list[[covariates[j]]] <- param_summary
  }
  
  par(mfrow = c(1, 1))
  summary_list
}


PosteriorDistributionSigmadelta <- function(mcmc_output, 
                                            true_sigma2 = NULL, 
                                            true_delta2 = NULL,
                                            trim_quantiles = c(0.01, 0.99)) {
  
  sigma2s <- as.numeric(unlist(mcmc_output$sigma2s))
  delta2s <- as.numeric(unlist(mcmc_output$delta2s))
  
  param_samples <- list(sigma2s = sigma2s, delta2s = delta2s)
  param_names <- c(expression(sigma^2), expression(delta^2))
  
  posterior_medians <- sapply(param_samples, median)
  credible_intervals <- lapply(param_samples, function(x)
    quantile(x, probs = c(0.025, 0.975)))
  
  par(mfrow = c(1, 2))
  
  for (j in 1:2) {
    samples <- param_samples[[j]]
    ci <- credible_intervals[[j]]
    true_val <- if (j == 1) true_sigma2 else true_delta2
    
    # Trim extreme outliers for visualization
    bounds <- quantile(samples, probs = trim_quantiles)
    samples_trimmed <- samples[samples >= bounds[1] & samples <= bounds[2]]
    
    # Define plotting range
    x_min <- min(samples_trimmed)* 0.99
    x_max <- max(samples_trimmed)* 1.01
    
    hist_info <- hist(samples_trimmed,
                      main = bquote("Posterior of " ~ .(param_names[[j]])),
                      xlab = as.expression(param_names[[j]]),
                      col = "gray",
                      border = "white",
                      breaks = 20,
                      freq = TRUE,
                      xlim = c(x_min, x_max))
    
    y_max <- max(hist_info$counts)
    
    # Credible interval
    rect(xleft = ci[1],
         xright = ci[2],
         ybottom = 0,
         ytop = y_max,
         col = rgb(0, 0, 0.2, alpha = 0.2),
         border = NA)
    
    # Posterior mean
    abline(v = posterior_medians[[j]], col = "blue", lwd = 2, lty = 2)
    
    # True value (if within trimmed range)
    if (!is.null(true_val) && true_val >= x_min && true_val <= x_max) {
      abline(v = true_val, col = "green", lwd = 2)
    }
  }
  
  par(mfrow = c(1, 1))
  
  # Return credible intervals nicely formatted
  credible_intervals <- do.call(rbind, credible_intervals)
  rownames(credible_intervals) <- c("sigma^2", "delta^2")
  print(credible_intervals)
  credible_intervals
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

TracePlots = function(sigma2_samples, delta2_samples, 
                      true_sigma2 = NULL, true_delta2 = NULL) {
  # Convert lists to numeric vectors
  sigma2_vec <- unlist(sigma2_samples)
  delta2_vec <- unlist(delta2_samples)
  
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # Two plots stacked
  
  # Trace plot for log(sigma²)
  plot(sigma2_vec, type = "l", col = "blue",
       ylab = expression(log(sigma^2)), xlab = "Iteration",
       main = expression("Traceplot of log(" ~ sigma^2 ~ ")"),
       log = "y")
  if (!is.null(true_sigma2)) {
    abline(h = true_sigma2, col = "darkgreen", lwd = 2, lty = 2)
  }
  
  # Trace plot for log(delta²)
  plot(delta2_vec, type = "l", col = "purple",
       ylab = expression(log(delta^2)), xlab = "Iteration",
       main = expression("Traceplot of log(" ~ delta^2 ~ ")"),
       log = "y")
  if (!is.null(true_delta2)) {
    abline(h = true_delta2, col = "darkgreen", lwd = 2, lty = 2)
  }
  
  par(mfrow = c(1, 1))  # Reset layout
}

ChangepointsLocationProbabilities <- function(allocation_vector_list, plot = TRUE) {
  # Flatten the changepoints, removing first and last elements from each vector
  all_changepoints <- unlist(lapply(allocation_vector_list, function(x) x[-c(1, length(x))]))
  
  # Define possible changepoint positions
  positions <- seq(min(unlist(allocation_vector_list)) + 1,
                   max(unlist(allocation_vector_list)) - 1)
  
  # Count how often each position is selected as a changepoint
  changepoint_counts <- table(factor(all_changepoints, levels = positions))
  
  # Compute changepoint probabilities
  changepoint_probs <- changepoint_counts / length(allocation_vector_list)
  
  # Plot if requested
  if (plot) {
    barplot(changepoint_probs,
            main = "Changepoint Location Probabilities",
            xlab = "Position",
            ylab = "Probability",
            col = "skyblue",
            border = NA,
            las = 2)
  }
  
  # Return as named list
  as.list(changepoint_probs)
}

OrderedChangepointsDistribution <- function(allocation_vector_list, n_changepoints = NULL) {
  # Extract internal changepoints (exclude first and last positions)
  internal_changepoints <- lapply(allocation_vector_list, function(x) x[-c(1, length(x))])
  
  # Count internal changepoints per allocation vector
  changepoint_counts <- sapply(internal_changepoints, length)
  
  # Filter by specific number of changepoints if provided
  if (!is.null(n_changepoints)) {
    internal_changepoints <- internal_changepoints[changepoint_counts == n_changepoints]
    if (length(internal_changepoints) == 0) {
      stop("No allocation vectors with the specified number of changepoints.")
    }
    changepoint_counts <- rep(n_changepoints, length(internal_changepoints))
  } else {
    changepoint_counts <- sapply(internal_changepoints, length)
  }
  
  # Determine full temporal range from all allocation vectors
  min_time <- min(unlist(allocation_vector_list))
  max_time <- max(unlist(allocation_vector_list))
  
  # Determine max number of changepoints to align row lengths
  max_changepoints <- max(sapply(internal_changepoints, length))
  
  # Pad with NA to get rectangular structure
  padded_changepoints <- lapply(internal_changepoints, function(x) {
    length(x) <- max_changepoints
    return(x)
  })
  
  # Build data frame: each row is one allocation vector
  changepoint_df <- as.data.frame(do.call(rbind, padded_changepoints))
  colnames(changepoint_df) <- paste0("Changepoint_", seq_len(max_changepoints))
  changepoint_df$n_changepoints <- changepoint_counts
  
  # Convert to long format for plotting
  changepoint_long <- reshape2::melt(changepoint_df, 
                                     id.vars = "n_changepoints", 
                                     variable.name = "Order", 
                                     value.name = "Position", 
                                     na.rm = TRUE)
  
  # Plot: horizontal boxplots with fixed x-axis limits
  p <- ggplot(changepoint_long, aes(x = Position, y = Order, fill = Order)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, outlier.alpha = 0.5) +
    scale_x_continuous(limits = c(min_time, max_time)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = "Ordered Changepoint Position Distributions",
      x = "Position (Time)",
      y = "Changepoint Order"
    )
  
  # Facet by number of changepoints if not filtering
  if (is.null(n_changepoints)) {
    p <- p + facet_wrap(~ n_changepoints, ncol = 1, labeller = label_both)
  }
  
  p
}




InclusionProbs <- function(models_list, covariates, decreasing = TRUE) {
  flat_models <- unlist(models_list)
  p <- length(models_list[[1]])
  inclusion_probs <- sapply(1:p, function(i) {
    mean(flat_models[seq(i,length(flat_models),by=p)])
  })
  # Assign names
  names(inclusion_probs) <- covariates
  
  # Sort in decreasing order
  inclusion_probs <- inclusion_probs[order(inclusion_probs, decreasing = decreasing)]
  
  round(inclusion_probs,3)
}

MultivariateNormalSample = function(N = 1, mu = 0, Sigma = 1){
  as.numeric(MASS::mvrnorm(N, mu = mu, Sigma = Sigma))
} 

MultivariateNormalDensity = function(X, mu = rep(0,length(X)), Sigma = 1, log = TRUE){
  if (length(Sigma) != 1){
    dmvnorm(x = X, mu = mu, Sigma = Sigma, log = log)
  }
  else{
    if (log) {
      as.numeric(- 0.5 * length(X) * log(2*pi) - 0.5 * length(X) * log(Sigma) - 0.5 / Sigma * sum((X - mu)^2))  
    } 
    else {
      as.numeric( 1 / (sqrt(2 * Sigma * pi)^(length(X))) * exp(- 1/(2*Sigma) * sum((X - mu)^2)))
    }  
  }
}

InverseGammaSample = function(N = 1, a = 0.01, b = 0.01){
  1 / rgamma(N, shape = a, rate = b)
} 

InverseGammaDensity = function(X, a = 0.01, b = 0.01, log = TRUE){
  if (log){
    a * log(b) - log(gamma(a)) - (a + 1) * log(X) - b/X
  } else {
    (b^a / gamma(a)) * X^(-a-1) * exp(-b/X)
  }
}

PoissonDensity <- function(X, lambda = 1, log = TRUE){
  if (log){
    X * log(lambda) - lambda - lgamma(X + 1)
  }
  else{
    lambda^X * exp(-lambda)/(factorial(X))
  }
}

NegativeBinomialDensity <- function(X, lambda = 1, r = 1, log = TRUE){
  if (log){
    LogChoose(n = X + r - 1, k = r - 1) + r * (log(r) - log(lambda + r)) + X * (log(lambda) - log(lambda + r))
  }
  else{
    choose(X + r - 1, r - 1) * (r/(lambda + r))^r * (lambda/(lambda + r))^X
  }
}


GeometricDensity <- function(X, rho = 0.01, log = TRUE){
  if (log){
    log(rho) + (X-1)*log(1-rho)
  }
  else{
    rho * (1-rho)^(X-1)
  }
}

CumulativeComplementary <- function(X, rho = 0.01, log = TRUE){
  ifelse(log, (X) * log(1-rho), (1-rho)^(X))
}

LogChoose <- function(n, k) {
  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
}

ModelUpdate <- function(k, model = rep(1, k), max_parents = length(model), intercept_move = FALSE) {
  model_star <- model
  idx_start <- 2 - intercept_move
  model_subset <- model[idx_start:k]
  model_sum <- sum(model)
  
  available_add <- which(model_subset == 0) + (1 - intercept_move)
  available_del <- which(model_subset == 1) + (1 - intercept_move)
  
  if (model_sum == 1 - intercept_move) {
    move_type <- 1  # ADD
    sampled <- sample(available_add, 1)
    model_star[sampled] <- 1
  } else if (model_sum == k - intercept_move) {
    move_type <- 2  # DELETE
    sampled <- sample(available_del, 1)
    model_star[sampled] <- 0
  } else {
    move_type <- ifelse(model_sum == max_parents, sample(2:3, 1), sample(1:3, 1))
    if (move_type == 1) {
      sampled <- if (length(available_add) == 1) available_add else sample(available_add, 1)
      model_star[sampled] <- 1
    } else if (move_type == 2) {
      sampled <- if (length(available_del) == 1) available_del else sample(available_del, 1)
      model_star[sampled] <- 0
    } else {  # EXCHANGE
      add_sample <- if (length(available_add) == 1) available_add else sample(available_add, 1)
      del_sample <- if (length(available_del) == 1) available_del else sample(available_del, 1)
      model_star[add_sample] <- 1
      model_star[del_sample] <- 0
    }
  }
  
  list(model = model_star, move = move_type)
}

BetasMeanVector <- function(Y, X, S,
                            model = rep(1, ncol(X)),
                            allocation_vector = c(1, max(S)),
                            init_betas = rep(0, ncol(X)),
                            sigma2 = 1, delta2 = 1, xi2 = 1, coupled = FALSE) {
  
  model_logical <- model == 1
  k <- sum(model_logical)
  H <- length(allocation_vector) - 1  # Number of segments
  betas_mean <- vector("list", H - 1)

  # Sanity check
  if (length(init_betas) != ncol(X)) {
    stop("init_betas must be of length equal to the number of columns in X (full parameter space).")
  }
  if (coupled){
    X_included <- X[, model_logical, drop = FALSE]
    init_betas_included <- init_betas[model_logical]
    
    # h = 1 → initial beta vector
    betas_mean[[1]] <- setNames(init_betas, colnames(X))  # already in full space
    
    # h = 2
    if (H >= 2) {
      X_curr <- as.matrix(X_included[S >= allocation_vector[2] & S < allocation_vector[3], , drop = FALSE])
      Y_curr <- Y[S >= allocation_vector[2] & S < allocation_vector[3]]
      
      beta_h2 <- solve((1 / delta2) * diag(k) + t(X_curr) %*% X_curr) %*%
        ((1 / delta2) * init_betas_included + t(X_curr) %*% Y_curr)
      
      beta_full <- rep(0, ncol(X))
      beta_full[model_logical] <- beta_h2[, 1]
      betas_mean[[2]] <- setNames(beta_full, colnames(X))
    }
    
    # h >= 3
    if (H >= 3) {
      for (h in 3:H) {
        X_curr <- as.matrix(X_included[S >= allocation_vector[h] & S < allocation_vector[h + 1], , drop = FALSE])
        Y_curr <- Y[S >= allocation_vector[h] & S < allocation_vector[h + 1]]
        
        beta_prev <- betas_mean[[h - 1]][model_logical]  # extract relevant subset
        
        beta_h <- solve((1 / xi2) * diag(k) + t(X_curr) %*% X_curr) %*%
          ((1 / xi2) * beta_prev + t(X_curr) %*% Y_curr)
        
        beta_full <- rep(0, ncol(X))
        beta_full[model_logical] <- beta_h[, 1]
        betas_mean[[h]] <- setNames(beta_full, colnames(X))
      }
    }
  }
  else{
    named_vector_template <- setNames(rep(0, ncol(X)), names(model))
    if (H == 1){
      betas_mean[[1]] <- named_vector_template
    }
    else{
      betas_mean <- replicate(H-1, named_vector_template, simplify = FALSE)
    }
  }
  return(betas_mean)
}

BinaryChangepoints <- function(changepoints, T = max(changepoints)) {
  internal_cp <- changepoints[changepoints > 1 & changepoints < T]
  binary_vec <- integer(T - 2)
  binary_vec[internal_cp - 1] <- 1
  binary_vec
}

ChangepointsList <- function(binary_vec) {
  internal_cps <- which(binary_vec == 1) + 1
  c(1, internal_cps, length(binary_vec) + 2)
}

ReallocationIndex <- function(x, i) {
  # Find left changepoint (exclusive)
  l <- i - 1
  while (l > 0 && x[l] == 0) l <- l - 1
  left_bound <- l + 1  # first non-changepoint to the right of previous CP
  
  # Find right changepoint (exclusive)
  r <- i + 1
  while (r <= length(x) && x[r] == 0) r <- r + 1
  right_bound <- r - 1  # last non-changepoint to the left of next CP
  
  candidates <- setdiff(left_bound:right_bound, i)
  if (length(candidates) == 0) {
    NaN
  }
  else if (length(candidates) == 1){
    candidates
  }  
  else{
    sample(candidates, 1)
  }
}


UpdateChangepoints <- function(allocation_vector, max_changepoints = (max(allocation_vector) - 2), n_segments = NaN) {
  allocation_vector_star <- BinaryChangepoints(allocation_vector)
  #print(allocation_vector_star)
  k <- length(allocation_vector_star)
  current_changepoints <- sum(allocation_vector_star)
  
  # Determine move type: 1 = BIRTH, 2 = DEATH, 3 = REALLOCATION
  if (!is.nan(n_segments)) {
    move_type <- 3  # Only reallocation allowed when n_segments is fixed
  } else if (current_changepoints == 0) {
    move_type <- 1  # Only BIRTH possible
  } else if (current_changepoints >= max_changepoints) {
    move_type <- sample(2:3, 1)  # Only DEATH or REALLOCATION
  } else {
    move_type <- sample(1:3, 1)  # Any move
  }
  
  # Perform the move
  if (move_type == 1) {  # BIRTH
    zero_indices <- which(allocation_vector_star == 0)
    if (length(zero_indices) > 0) {
      sampled_changepoint <- if (current_changepoints == (k - 1)) {
        zero_indices
      } else {
        sample(zero_indices, 1)
      }
      allocation_vector_star[sampled_changepoint] <- 1
    }
  } else if (move_type == 2) {  # DEATH
    one_indices <- which(allocation_vector_star == 1)
    if (length(one_indices) > 0) {
      sampled_changepoint <- if (current_changepoints == 1) {
        one_indices
      } else {
        sample(one_indices, 1)
      }
      allocation_vector_star[sampled_changepoint] <- 0
    }
  } else if (move_type == 3) {  # REALLOCATION
    one_indices <- which(allocation_vector_star == 1)
    if (length(one_indices) > 0) {
      sampled_changepoint <- if (length(one_indices) == 1) {
        one_indices
      } else {
        sample(one_indices, 1)
      }
      sampled_changepoint_2 <- ReallocationIndex(allocation_vector_star, sampled_changepoint)
      if (!is.nan(sampled_changepoint_2)) {
        allocation_vector_star[sampled_changepoint] <- 0
        allocation_vector_star[sampled_changepoint_2] <- 1
      }
    }
  }
  
  allocation_vector <- ChangepointsList(allocation_vector_star)
  list(allocation_vector = allocation_vector, move = move_type)
}

BetasFCDSample <- function(N = 1, Y, X, segment_indexes, sigma2, delta2, mu_betas = 0, sigma2_in_betas_prior = FALSE) {
  # Loop through segments based on allocation_vector and S
  betas_list <- lapply(1:max(segment_indexes), function(h) {
    X_current <- X[segment_indexes == h, , drop = FALSE]
    Y_current <- Y[segment_indexes == h]
    # Compute the posterior for each segment
    #V_current <- if (sigma2_in_betas_prior) {
    #  solve(t(X_current) %*% X_current + 1/delta2 * diag(ncol(X_current)))
    #} else {
    #  solve(1/sigma2 * t(X_current) %*% X_current + 1/delta2 * diag(ncol(X_current)))
    #}
    
    V_current <- tryCatch({
      A <- if (sigma2_in_betas_prior) {
        t(X_current) %*% X_current + 1 / delta2 * diag(ncol(X_current))
      } else {
        1 / sigma2 * t(X_current) %*% X_current + 1 / delta2 * diag(ncol(X_current))
      }
      if (rcond(A) < 1e-10) {
        #message("Matrix is near-singular, using pseudo-inverse.")
        ginv(A)  # fallback to Moore-Penrose pseudo-inverse
      } else {
        solve(A)
      }
    }, error = function(e) {
      #message("solve() failed, using pseudo-inverse.")
      ginv(A)
    })
    
    mu_betas_h <- if (sigma2_in_betas_prior) {
      V_current %*% (t(X_current) %*% Y_current + 1/delta2 * mu_betas[[h]])
    } else {
      V_current %*% (1/sigma2 * t(X_current) %*% Y_current + 1/delta2 * mu_betas[[h]])
    }
    # Sample betas for this segment
    MultivariateNormalSample(N = N, mu = mu_betas_h, Sigma = if (sigma2_in_betas_prior) sigma2 * V_current else V_current)
  })
  
  betas_list 
}

BetasFCDDensity = function(Y, X, betas = rep(0, ncol(X)), sigma2 = 1, delta2 = 1, mu_betas = 0, sigma2_in_betas_prior = FALSE, log = TRUE, segment_indexes = rep(1, length(Y))){
  # Loop through segments based on allocation_vector and S
  betas_list <- sapply(1:max(segment_indexes), function(h) {
    X_current <- X[segment_indexes == h, , drop = FALSE]
    Y_current <- Y[segment_indexes == h]
    betas_current <- betas[[h]]
    
    # Compute the posterior for each segment
    V_current <- if (sigma2_in_betas_prior) {
      solve(t(X_current) %*% X_current + 1/delta2 * diag(ncol(X_current)))
    } else {
      solve(1/sigma2 * t(X_current) %*% X_current + 1/delta2 * diag(ncol(X_current)))
    }
    
    mu_betas_h <- if (sigma2_in_betas_prior) {
      V_current %*% (t(X_current) %*% Y_current + 1/delta2 * mu_betas[[h]])
    } else {
      V_current %*% (1/sigma2 * t(X_current) %*% Y_current + 1/delta2 * mu_betas[[h]])
    }

    # Sample betas for this segment
    MultivariateNormalDensity(X = betas_current, mu = t(mu_betas_h), Sigma = if (sigma2_in_betas_prior) sigma2 * V_current else V_current, log = log)
  })
  sum(betas_list)
}

Sigma2FCDSample <- function(N = 1, Y, X, betas, segment_indexes,
                            delta2 = 1, mu_betas = 0, a = 0.01, b = 0.01,
                            sigma2_in_betas_prior = FALSE, collapsed = FALSE) {
  
  # No need for 'model' logic, so we directly use the full X and betas
  k <- ncol(X)  # Number of columns in X (features)
  
  total <- 0  # Initialize the total sum
  
  for (h in 1:max(segment_indexes)) {
    X_current <- X[segment_indexes == h, , drop = FALSE]  # Extract data for the group
    Y_current <- Y[segment_indexes == h]  # Corresponding response values for the group
    betas_current <- betas[[h]]  # Coefficients for this group
    
    if (!sigma2_in_betas_prior) {
      res <- Y_current - X_current %*% betas_current  # Residuals for the group
      total <- total + sum(res^2)  # Sum of squared residuals
    } 
    else if (!collapsed) {
      res <- Y_current - X_current %*% betas_current  # Residuals for the group
      res_beta <- betas_current - mu_betas[[h]]  # Difference between coefficients and prior mean
      total <- total + sum(res^2) + sum(res_beta^2) / delta2  # Add to total
    } 
    else {
      D <- diag(1 / delta2, k)  # Prior precision for coefficients
      M <- t(X_current) %*% Y_current + D %*% mu_betas[[h]]  # Posterior mean calculation
      V <- solve(t(X_current) %*% X_current + D, M)  # Posterior covariance
      total <- total + as.numeric(t(Y_current) %*% Y_current - t(M) %*% V)  # Update total
    }
  }
  
  # Sampling from the inverse gamma distribution
  InverseGammaSample(N = N, a = a + 0.5 * length(Y), b = b + 0.5 * total)
}



Delta2FCDSample <- function(N = 1, betas_list, segment_indexes,
                            sigma2 = 1, alpha = 0.01, beta = 0.01,
                            mu_betas = 0, sigma2_in_betas_prior = FALSE) {
  
  k <- length(mu_betas[[1]])
  total <- 0
  
  for (h in 1:max(segment_indexes)) {
    res <- betas_list[[h]] - mu_betas[[h]]
    total <- total + if (!sigma2_in_betas_prior) sum(res^2) else sum(res^2) / sigma2
  }

  InverseGammaSample(N = N, a = alpha + 0.5 * k, b = beta + 0.5 * total)
}

Likelihood <- function(Y, X, betas = rep(0, ncol(X)), sigma2 = 1, segment_indexes = rep(1, length(Y)), log = TRUE){
  sum(sapply(1:max(segment_indexes), function(i) {MultivariateNormalDensity(X = Y[segment_indexes==i], mu = X[segment_indexes==i, ,drop = FALSE] %*% betas[[i]], Sigma = sigma2, log = log)}))
}

BetaPrior <- function(betas, betas_0 = rep(0, length(betas[[1]])), segment_indexes = rep(1, length(Y)), sigma2 = 1, delta2 = 1, sigma2_in_betas_prior = FALSE, log = TRUE){
  BetaNoise = ifelse(sigma2_in_betas_prior, sigma2*delta2, delta2)
  sum(sapply(1:max(segment_indexes), function(i) {MultivariateNormalDensity(betas[[i]], mu = betas_0[[i]], Sigma = BetaNoise, log = log)}))
}

#LogMarginal = function(Y, X, delta2 = 1, mu_beta = rep(0, ncol(X)), a = 0.01, b = 0.01, segment_indexes= rep(1, length(Y))){
#  k <- ncol(X)
#  N <- length(Y)
#  D <- sum(sapply(1:max(segment_indexes), function(h) {
#    X_current <- X[segment_indexes == h, , drop = FALSE]
#    as.numeric(determinant(diag(k) + delta2 * t(X_current) %*% X_current, logarithm = TRUE)$modulus)
#    }))
#  E <- sum(sapply(1:max(segment_indexes), function(h) {
#    X_current <- X[segment_indexes == h, , drop = FALSE]  # Extract data for the group
#    Y_current <- Y[segment_indexes == h]  # Corresponding response values for the group
#    mu_beta_current <- mu_beta[[h]]
#    res <- (Y_current - X_current %*% mu_beta_current)
#    M <- (diag(length(Y_current)) - delta2 * X_current %*% solve(diag(k) + delta2 * t(X_current) %*% X_current) %*% t(X_current))
#    t(res) %*% M %*% res
#  }))
#  as.numeric(lgamma(N/2 + a) - lgamma(a) - N/2 * log(pi) + a * log(2*b) - D/2 - (N/2 + a) * log(2*b + E))
#}
LogMarginal <- function(Y, X, delta2 = 1, mu_beta = rep(0, ncol(X)),
                        a = 0.01, b = 0.01, segment_indexes = rep(1, length(Y))) {
  library(MASS)
  k <- ncol(X)
  N <- length(Y)
  
  D <- sum(sapply(1:max(segment_indexes), function(h) {
    tryCatch(
      as.numeric(determinant(diag(k) + delta2 * t(X[segment_indexes == h, , drop = FALSE]) %*% 
                               X[segment_indexes == h, , drop = FALSE], logarithm = TRUE)$modulus),
      error = function(e) {
        message("Falling back to log-det with ridge for segment ", h)
        log(det(diag(k) + delta2 * t(X[segment_indexes == h, , drop = FALSE]) %*% 
                  X[segment_indexes == h, , drop = FALSE] + 1e-6 * diag(k)))
      }
    )
  }))
  
  ridge <- 1e-6
  E <- sum(sapply(1:max(segment_indexes), function(h) {
    Xh <- X[segment_indexes == h, , drop = FALSE]
    Yh <- Y[segment_indexes == h]
    Ah <- diag(k) + delta2 * t(Xh) %*% Xh + ridge * diag(k)
    inv_Ah <- tryCatch(
      chol2inv(chol(Ah)),
      error = function(e) {
        message(paste("chol2inv failed in segment", h, "- falling back to pseudo-inverse"))
        ginv(Ah)
      }
    )
    Mh <- diag(nrow(Xh)) - delta2 * Xh %*% inv_Ah %*% t(Xh)
    res <- Yh - Xh %*% mu_beta[[h]]
    as.numeric(t(res) %*% Mh %*% res)
  }))
  E <- max(E,0)
  as.numeric(lgamma(N / 2 + a) - lgamma(a) - N / 2 * log(pi) +
               a * log(2 * b) - D / 2 - (N / 2 + a) * log(2 * b + E))
}



ModelPrior <- function(model, type = 'Uniform', max_parents = length(model) - 1 + intercept_move, lambda = 1, r = 1, intercept_move = FALSE, log = TRUE){
  if (sum(model) > max_parents){
    ifelse(log, - Inf, 0)
  }
  else{
    if (type == 'Uniform'){ # Uniform over the models with |parent_set| <= max_parents
      ifelse(log, 0, 1)
    }
    else if (type == 'Uniform_card'){ # Uniform over the models with the same cardinality and over the cardinalities
      ifelse(log, - log(max_parents + 1) - LogChoose(max_parents, (sum(model)) - 1 + intercept_move), 1/((max_parents + 1) * choose(max_parents,(sum(model)) - 1 + intercept_move)))
    }
    else if (type == 'Poisson'){ # Poisson distribution of the parent set's cardinality 
      PoissonDensity(X = (sum(model)-1+ intercept_move), lambda = lambda, log = log)
    }
    else if (type == 'Neg-Bin'){ # Negative Binomial distribution of the parent set's cardinality 
      NegativeBinomialDensity(X = (sum(model)-1+ intercept_move), lambda = lambda, r = r, log = log)
    }
  }
}

SegmentPrior <- function(allocation_vector, type_segment_prior = 'Uniform_card', rho = 0.01, iota = 1, log = TRUE, max_changepoints = max(allocation_vector)){
  if ((length(allocation_vector) - 2) > max_changepoints){
    ifelse(log, - Inf, 0)
  }
  else{
    if (type_segment_prior == 'Uniform_card'){
      ifelse(log,0,1)
    }
    else if (type_segment_prior == 'Geometric'){
      distances <- diff(allocation_vector)
      if (log){
        ifelse(length(distances) == 1,GeometricDensity(X = distances, rho = rho, log = TRUE),sum(sapply(1:(length(distances)-1), function(i) {GeometricDensity(X = distances[i], rho = rho, log = TRUE)}))) + CumulativeComplementary(X = distances[length(distances)], rho = rho, log = TRUE)
      }
      else{
        ifelse(length(distances) == 1,GeometricDensity(X = distances, rho = rho, log = FALSE),prod(sapply(1:(length(distances)-1), function(i) {GeometricDensity(X = distances[i], rho = rho, log = FALSE)}))) * CumulativeComplementary(X = distances[length(distances)], rho = rho, log = FALSE)
      }
    }
    else if (type_segment_prior == 'Poisson_card'){
      PoissonDensity(length(allocation_vector) - 2, lambda = iota, log = log)
    }
  }
}

PosteriorRatio = function(Y, X, X_star = X, model = rep(1, ncol(X)), model_star = model, betas = rep(0, ncol(X)), betas_star = betas, sigma2 = 1, sigma2_star = sigma2, delta2 = 1, delta2_star = delta2, betas_0 = rep(0, ncol(X)), betas_0_star = betas_0, a = 0.01, b = 0.01, alpha = 0.01, beta = 0.01, lambda = 1, r = 1, rho = 0.01, iota = 1, type_model_prior = 'Uniform', type_segment_prior = 'Uniform_card', intercept_move = intercept_move, max_parents = length(model) - 1 + intercept_move, max_changepoints = 0, likelihood = TRUE, marginal = FALSE, model_prior = FALSE, segment_prior = FALSE, beta_prior = TRUE, sigma_prior = FALSE, delta_prior = FALSE, log = TRUE, sigma2_in_betas_prior = FALSE, segment_indexes = rep(1, length(Y)), segment_indexes_star = segment_indexes, allocation_vector = c(1, max(1)), allocation_vector_star = allocation_vector){
  if (log){
    MarginalRate = ifelse(marginal,LogMarginal(Y = Y, X = X_star, delta2 = delta2, mu_beta = betas_0_star, a = a, b = b, segment_indexes = segment_indexes_star) - LogMarginal(Y = Y, X = X, delta2 = delta2, mu_beta = betas_0, a = a, b = b, segment_indexes = segment_indexes),0)
    LikelihoodRate = ifelse(likelihood,Likelihood(Y = Y, X = X_star, betas = betas_star, sigma2 = sigma2_star, log = TRUE, segment_indexes = segment_indexes_star) - Likelihood(Y = Y, X = X, betas = betas, sigma2 = sigma2, log = TRUE, segment_indexes = segment_indexes),0)
    BetaPriorRate = ifelse(beta_prior, BetaPrior(betas = betas_star, betas_0 = betas_0_star, sigma2 = sigma2_star, delta2 = delta2_star, sigma2_in_betas_prior = sigma2_in_betas_prior, log = TRUE, segment_indexes = segment_indexes_star) - BetaPrior(betas = betas, betas_0 = betas_0, sigma2 = sigma2, delta2 = delta2, sigma2_in_betas_prior = sigma2_in_betas_prior, log = TRUE, segment_indexes = segment_indexes), 0)
    SigmaPriorRate = ifelse(sigma_prior, InverseGammaDensity(sigma2_star, a = a, b = b, log = TRUE) - InverseGammaDensity(sigma2, a = a, b = b, log = TRUE), 0)
    DeltaPriorRate = ifelse(delta_prior, InverseGammaDensity(delta2_star, a = alpha, b = beta, log = TRUE) - InverseGammaDensity(delta2, a = alpha, b = beta, log = TRUE), 0)
    ModelPriorRate = ifelse(model_prior, ModelPrior(model_star, type = type_model_prior, lambda = lambda, r = r, max_parents = max_parents, intercept_move = intercept_move, log = TRUE) - ModelPrior(model, type = type_model_prior, lambda = lambda, r = r, max_parents = max_parents, intercept_move = intercept_move, log = TRUE), 0)
    SegmentPriorRate = ifelse(segment_prior, SegmentPrior(allocation_vector = allocation_vector_star, type_segment_prior = type_segment_prior, rho = rho, iota = iota, max_changepoints = max_changepoints, log = TRUE) - SegmentPrior(allocation_vector = allocation_vector, type_segment_prior = type_segment_prior, rho = rho, iota = iota, max_changepoints = max_changepoints, log = TRUE),0)
    PosteriorRate = MarginalRate + LikelihoodRate + BetaPriorRate + SigmaPriorRate + DeltaPriorRate + ModelPriorRate + SegmentPriorRate
  }
  else{
    MarginalRate = ifelse(marginal,Marginal(Y = Y, X = X_star, delta2 = delta2, mu_beta = betas_0_star, a = a, b = b) / Marginal(Y = Y, X = X, delta2 = delta2, mu_beta = betas_0, a = a, b = b),1)
    LikelihoodRate = ifelse(likelihood,MultivariateNormalDensity(Y, mu = X_star %*% betas_star, Sigma = sigma2_star, log = FALSE) / MultivariateNormalDensity(Y, mu = X %*% betas, Sigma = sigma2, log = FALSE),1)
    BetaPriorRate = ifelse(beta_prior, MultivariateNormalDensity(betas_star, mu = betas_0, Sigma = BetaNoise_star, log = FALSE) / MultivariateNormalDensity(betas, mu = betas_0, Sigma = BetaNoise, log = FALSE), 1)  
    SigmaPriorRate = ifelse(sigma_prior, InverseGammaDensity(sigma2_star, a = a, b = b, log = FALSE) / InverseGammaDensity(sigma2, a = a, b = b, log = FALSE), 1)
    DeltaPriorRate = ifelse(delta_prior, InverseGammaDensity(delta2_star, a = alpha, b = beta, log = FALSE) / InverseGammaDensity(delta2, a = alpha, b = beta, log = FALSE), 1)
    ModelPriorRate = ifelse(model_prior, ModelPrior(model_star, type = type_model_prior, lambda = lambda, r = r, max_parents = max_parents, intercept_move = intercept_move, log = FALSE) / ModelPrior(model, type = type_model_prior, lambda = lambda, r = r, max_parents = max_parents, intercept_move = intercept_move, log = FALSE), 1)
    SegmentPriorRate = ifelse(segment_prior, SegmentPrior(allocation_vector = allocation_vector_star, type_segment_prior = type_segment_prior, rho = rho, iota = iota, log = FALSE) / SegmentPrior(allocation_vector = allocation_vector, type_segment_prior = type_segment_prior, rho = rho, iota = iota, log = FALSE),1)
    PosteriorRate = LikelihoodRate * BetaPriorRate * SigmaPriorRate * DeltaPriorRate * ModelPriorRate * SegmentsPriorRate
  }
  PosteriorRate
}

ProposalRatio <- function(Y, X, X_star = X, model = rep(1, ncol(X)), betas = rep(0, ncol(X)), betas_star = betas, sigma2 = 1, delta2 = 1, betas_0 = rep(0, ncol(X)), betas_0_star = betas_0, type = 'Gibbs', move_type = 1, max_parents = length(model) - 1 + intercept_move, intercept_move = FALSE, log = TRUE, segment_indexes = rep(1, length(Y)), segment_indexes_star = segment_indexes, allocation_vector = c(1, 1)){
  p <- sum(model)
  k <- length(model)
  TT <- max(allocation_vector)
  tau_card <- length(allocation_vector) - 2
  if (type %in% c('MH','Gibbs')){
    if (log){
      0
    }
    else{
      1
    }
  }
  else if (type %in% c('RJ-FCD')){
    if (move_type == 1){ # addition
      cov_proposal <- ifelse(log, log(k - p) - log(p - intercept_move), (k - p) / (p - intercept_move))
    }
    else if (move_type == 2){ # deletion
      cov_proposal <- ifelse(log, log(p - 1 + intercept_move) - log(k - p + 1), (p - 1 + intercept_move) / (k - p + 1))
    }
    else{ # exchange
      cov_proposal <- 0
    }
    # Log-posterior ratio (likelihood + prior + proposal ratio)
    (- BetasFCDDensity(Y = Y, X = X_star, betas = betas_star, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_0_star, sigma2_in_betas_prior = FALSE, log = log, segment_indexes = segment_indexes_star) +
        BetasFCDDensity(Y = Y, X = X, betas = betas, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_0, sigma2_in_betas_prior = FALSE, log = log, segment_indexes = segment_indexes) + cov_proposal)
  }
  else if (type %in% c('RJ-Prior','RJ-Marginal')){
    if (move_type == 1){ # addition
      ifelse(log, log(k - p) - log(p - intercept_move), (k - p) / (p - intercept_move))
    }
    else if (move_type == 2){ # deletion
      ifelse(log, log(p - 1 + intercept_move) - log(k - p + 1), (p - 1 + intercept_move) / (k - p + 1))
    }
    else{ # reversal
      ifelse(log,0,1)
    }
  }
  else if (type %in% c('ChangePoint')){
    if (move_type == 1){
      ifelse(log, log(TT - tau_card) - log(tau_card - intercept_move), (TT - tau_card) / (tau_card - intercept_move))
    }
    else if (move_type == 2){
      ifelse(log, log(TT - tau_card) - log(tau_card - intercept_move), (TT - tau_card) / (tau_card - intercept_move))
    }
    else{ # reallocation
      ifelse(log,0,1)
    }
  }
}

AcceptanceRatio = function(Y, X, X_star = X, model = rep(1, ncol(X)), model_star = model, MCMC_type = 'Gibbs', betas = rep(0, ncol(X)), betas_star = rep(0, ncol(X_star)), sigma2 = 1, sigma2_star = 1, delta2 = 1, delta2_star = 1, betas_0 = rep(0, ncol(X)), betas_0_star = betas_0, a = 0.01, b = 0.01, alpha = 0.01, beta = 0.01, lambda = 1, r = 1, rho = 0.01, iota = 1, likelihood = TRUE, marginal = FALSE, model_prior = FALSE, segment_prior = FALSE, beta_prior = TRUE, sigma_prior = FALSE, delta_prior = FALSE, log = TRUE, sigma2_in_betas_prior = FALSE, move_type = 1, type_model_prior = 'Uniform', type_segment_prior = 'Uniform_card',intercept_move = FALSE, max_parents = length(model) - 1 + intercept_move, max_changepoints = 0, segment_indexes = rep(1, length(Y)), segment_indexes_star = segment_indexes, allocation_vector = c(1, max(1)), allocation_vector_star = allocation_vector){
  if (log){
    exp(PosteriorRatio(Y = Y, X = X, X_star = X_star, model = model, model_star = model_star, betas = betas, betas_star = betas_star, sigma2 = sigma2, sigma2_star = sigma2_star, delta2 = delta2, delta2_star = delta2_star, betas_0 = betas_0, betas_0_star = betas_0_star, a = a, b = b, alpha = alpha, beta = beta, lambda = lambda, r = r, likelihood = likelihood, marginal = marginal, model_prior = model_prior, beta_prior = beta_prior, sigma_prior = sigma_prior, delta_prior = delta_prior, log = TRUE, sigma2_in_betas_prior = sigma2_in_betas_prior, type_model_prior = type_model_prior, intercept_move = intercept_move, max_parents = max_parents, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes_star, allocation_vector = allocation_vector, allocation_vector_star = allocation_vector_star, segment_prior = segment_prior, rho = rho, iota = iota, type_segment_prior = type_segment_prior, max_changepoints = max_changepoints) + ProposalRatio(Y = Y, X = X, X_star = X_star, model = model, betas = betas, betas_star = betas_star, sigma2 = sigma2, delta2 = delta2, betas_0 = betas_0, betas_0_star = betas_0_star, type = MCMC_type, move_type = move_type, intercept_move = intercept_move, max_parents = max_parents,log = TRUE, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes_star, allocation_vector = allocation_vector))
  }  
  else{
    PosteriorRatio(Y = Y, X = X, X_star = X_star, model = model, model_star = model_star, betas = betas, betas_star = betas_star, sigma2 = sigma2, sigma2_star = sigma2_star, delta2 = delta2, delta2_star = delta2_star, betas_0 = betas_0, betas_0_star = betas_0_star, a = a, b = b, alpha = alpha, beta = beta, lambda = lambda, r = r, likelihood = likelihood, marginal = marginal, model_prior = model_prior, beta_prior = beta_prior, sigma_prior = sigma_prior, delta_prior = delta_prior, log = FALSE, sigma2_in_betas_prior = sigma2_in_betas_prior, type_model_prior = type_model_prior, intercept_move = intercept_move, max_parents = max_parents, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes_star, allocation_vector = allocation_vector, allocation_vector_star = allocation_vector_star, segment_prior = segment_prior, rho = rho, iota = iota, type_segment_prior = type_segment_prior, max_changepoints = max_changepoints) * ProposalRatio(Y = Y, X = X, X_star = X_star, model = model, betas = betas, betas_star = betas_star, sigma2 = sigma2, delta2 = delta2, betas_0 = betas_0, betas_0_star = betas_0_star, type = MCMC_type, move_type = move_type, intercept_move = intercept_move, max_parents = max_parents, log = FALSE, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes_star,allocation_vector = allocation_vector)
  }
}

MH_move = function(Y, X, S, model = rep(1,ncol(X)), betas = rep(0,ncol(X)), sigma2 = 1, delta2 = 1, alpha = 0.01, beta = 0.01, a = 0.01, b = 0.01, betas_mean = rep(0,ncol(X)), sigma2_in_betas_prior = FALSE, v = 0.1, log = TRUE, allocation_vector = c(1, max(S)), coupling = FALSE){
  betas_included <- lapply(if (is.list(betas)) betas else list(betas), function(b) b[model == 1])
  betas_mean <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  X = as.matrix(X[,model == 1])
  segment_indexes <- if(length(allocation_vector) > 2){
    cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    S
  }
  # Betas update
  for (j in 1:max(segment_indexes)){
    for (i in 1:sum(model)) {
      betas_star_included <- betas_included
      betas_star_included[[j]][i] <- betas_included[[j]][i] + runif(1, min = -v, max = v)
      
      A = AcceptanceRatio(Y = Y, X = X, X_star = X, MCMC_type = 'MH', betas = betas_included, betas_star = betas_star_included, delta2 = delta2, delta2_star = delta2, sigma2 = sigma2, sigma2_star = sigma2, betas_0 = betas_mean, betas_0_star = betas_mean, log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes)
      #cat(0.5 * (1 / sigma2) * (sum((Y - as.matrix(X) %*% betas) ^ 2) - sum((Y - as.matrix(X) %*% betas_star) ^ 2)) + 0.5 * 1 / delta2 * (betas[i] ^2 - betas_star[i] ^ 2), '\n')
      if (min(A,1) >= runif(1)) {
        betas_included[[j]] <- betas_star_included[[j]]
      }
    }
  }
  
  # Sigma update
  sigma2_star = abs(sigma2 + runif(1, min = -v, max = v))
  
  A = AcceptanceRatio(Y = Y, X = X, X_star = X, MCMC_type = 'MH', betas = betas_included, betas_star = betas_included, delta2 = delta2, delta2_star = delta2, sigma2 = sigma2, sigma2_star = sigma2_star, betas_0 = betas_mean, betas_0_star = betas_mean, a = a, b = b, beta_prior = sigma2_in_betas_prior, sigma_prior = TRUE, log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes)
  
  # Probability of Parameter Update: sigma^(t) <- sigma_star
  if (min(A,1) >= runif(1)) {
    sigma2 <- sigma2_star
  }
  
  # delta update
  delta2_star = abs(delta2 + runif(1, min = -v, max = v))
  
  A = AcceptanceRatio(Y = Y, X = X, X_star = X, MCMC_type = 'MH', betas = betas_included, betas_star = betas_included, delta2 = delta2, delta2_star = delta2_star, sigma2 = sigma2, sigma2_star = sigma2, betas_0 = betas_mean, betas_0_star = betas_mean, alpha = alpha, beta = beta, likelihood = FALSE, delta_prior = TRUE, log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes)
  
  # Probability of Parameter Update: delta^(t) <- delta_star
  if (min(A,1) >= runif(1)) {
    delta2 <- delta2_star
  }
  
  #betas <- replace(numeric(length(betas)), model == 1, betas_included)
  if (is.list(betas_included)) {
    betas <- lapply(betas_included, function(b_in) replace(numeric(length(model)), model == 1, b_in))
  } else {
    betas <- replace(numeric(length(model)), model == 1, betas_included)
  }
  return(list(betas = betas, sigma2 = sigma2, delta2 = delta2))
}

Gibbs_move = function(Y, X, S = rep(1, length(Y)), model = rep(1,ncol(X)), betas = rep(0,ncol(X)), sigma2 = 1, delta2 = 1, alpha = 0.01, beta = 0.01, a = 0.01, b = 0.01, betas_mean = rep(0,ncol(X)), sigma2_in_betas_prior = FALSE, collapsed = FALSE, allocation_vector = c(1, max(S)), coupling = FALSE){
  #betas_included <- lapply(betas, function(b) b[model == 1])
  betas_included <- lapply(if (is.list(betas)) betas else list(betas), function(b) b[model == 1])
  X <- as.matrix(X[,model == 1])
  betas_mean  <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  segment_indexes <- if(length(allocation_vector) > 2){
    cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    rep(1, length(Y))
  }
  # Sample sigma2 | betas, Y or sigma2 | Y
  sigma2 <- Sigma2FCDSample(N = 1, Y = Y, X = X, betas = betas_included, delta2 = delta2, a = a, b = b, mu_betas = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, collapsed = collapsed, segment_indexes = segment_indexes)
  # Sample betas | sigma2, delta2, Y
  betas_included <- BetasFCDSample(N = 1, Y = Y, X = X, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes)
  # Sample delta2 | betas
  delta2 <- Delta2FCDSample(N = 1, betas = betas_included, sigma2 = sigma2, alpha = alpha, beta = beta, mu_betas = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes)

  if (coupling){
    xi2 <- 1 # Add also xi2 update for delta2
  }
  
  #betas <- replace(numeric(length(betas)), model == 1, betas_included)
  if (is.list(betas_included)) {
    betas <- lapply(betas_included, function(b_in) replace(numeric(length(model)), model == 1, b_in))
  } else {
    betas <- replace(numeric(length(model)), model == 1, betas_included)
  }
  
  
  list(betas = betas, sigma2 = sigma2, delta2 = delta2)
}

RJ_MH_move = function(Y, X, S = rep(1, length(Y)), model = rep(1,ncol(X)), betas = rep(0,ncol(X)), sigma2 = 1, delta2 = 1, alpha = 0.01, beta = 0.01, a = 0.01, b = 0.01, r = 1, lambda = 1, betas_mean = rep(0,ncol(X)), sigma2_in_betas_prior = FALSE, log = TRUE, intercept_move = FALSE, type_model_prior = 'Uniform', max_parents = ncol(X), allocation_vector = c(1, max(S)), coupling = FALSE){
  k = ncol(X)
  model_star = model
  betas_star = betas
  betas_mean_star = betas_mean
  X_star = X
  X_included = as.matrix(X[,model == 1])
  betas_mean_included <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  betas_included = lapply(if (is.list(betas)) betas else list(betas), function(b) b[model == 1])
  segment_indexes <- if(length(allocation_vector) > 2){
    cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    rep(1, length(Y))
  }
  # Sample a random covariate X_j from the covariate vector
  sampled_cov <- ifelse(sum(model)==max_parents, (1 - intercept_move) + sample(which(model[(2 - intercept_move):k] == 1), 1), sample((2 - intercept_move):k, 1))
  move_type = ifelse(model[sampled_cov] == 1, 2, 1)
  # Covariate X_j addition/deletion from the model at time t-1
  model_star[sampled_cov] = 1 - model[sampled_cov]
  
  # Random sample of beta_j in case of covariate addition from a Normal distribution with delta^2 noise
  for (i in 1:max(segment_indexes)){
    betas_star[[i]][sampled_cov] = model_star[sampled_cov] * MultivariateNormalSample(1, mu = betas_mean[[i]][sampled_cov], Sigma = delta2)
  }
  X_star = as.matrix(X_star[,model_star == 1])
  betas_star_included = lapply(if (is.list(betas_star)) betas_star else list(betas_star), function(b) b[model_star == 1])
  betas_mean_star <- lapply(if (is.list(betas_mean_star)) betas_mean_star else list(betas_mean_star), function(b) b[model_star == 1])
  # Log Posterior Computation
  A = AcceptanceRatio(Y = Y, X = X_included, X_star = X_star, model = model, model_star = model_star, MCMC_type = 'RJ-Prior', betas = betas_included, betas_star = betas_star_included, delta2 = delta2, delta2_star = delta2, sigma2 = sigma2, sigma2_star = sigma2, betas_0 = betas_mean_included, betas_0_star = betas_mean_star, a = a, b = b, r = r, lambda = lambda, beta_prior = FALSE, model_prior = ifelse(type_model_prior=='Uniform', FALSE, TRUE), log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, intercept_move = intercept_move, type_model_prior = type_model_prior, max_parents = max_parents, move_type = move_type, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes)
  
  # Probability of Model Update: beta_j^(t) <- beta_j^(star) (X_j addition or deletion)
  if (min(A,1) >= runif(1, min = 0, max = 1)) {
    model <- model_star
    betas <- betas_star
  }
  
  list(betas = betas, model = model)
}

RJ_MH_FCD_move = function(Y, X, S = rep(1, length(Y)), model = rep(1,ncol(X)), betas = rep(0,ncol(X)), sigma2 = 1, delta2 = 1, alpha = 0.01, beta = 0.01, a = 0.01, b = 0.01, r = 1, lambda = 1, betas_mean = rep(0,ncol(X)), sigma2_in_betas_prior = FALSE, log = TRUE, intercept_move = FALSE, type_model_prior = 'Uniform', max_parents = ncol(X), allocation_vector = c(1, max(S)), coupling = FALSE){
  # Init parameters
  X_star <- X
  betas_mean_star <- betas_mean
  betas_mean <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  betas_included <- lapply(if (is.list(betas)) betas else list(betas), function(b) b[model == 1])
  
  segment_indexes <- if(length(allocation_vector) > 2){
    cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    rep(1, length(Y))
  }
  
  # Random covariate sampling (respecting intercept constraints)
  updated_model <- ModelUpdate(k = ncol(X), model = model, max_parents = max_parents, intercept_move = intercept_move)
  model_star <- updated_model$model
  move_type <- updated_model$move
  X <- as.matrix(X[,model==1])
  
  # Design matrices for current and proposed models
  X_star <- as.matrix(X_star[, model_star == 1])
  betas_mean_star <- lapply(if (is.list(betas_mean_star)) betas_mean_star else list(betas_mean_star), function(b) b[model_star == 1])
  
  # Prior precision matrix and FCD for betas_star
  betas_star_included <- BetasFCDSample(N = 1, Y = Y, X = X_star, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_mean_star, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes)
  # Acceptance Ratio (likelihood + prior + proposal ratio)
  A <- AcceptanceRatio(Y = Y, X = X, X_star = X_star, model = model, model_star = model_star, MCMC_type = 'RJ-FCD', betas = betas_included, betas_star = betas_star_included, delta2 = delta2, delta2_star = delta2, sigma2 = sigma2, sigma2_star = sigma2, betas_0 = betas_mean, betas_0_star = betas_mean_star, model_prior = ifelse(type_model_prior=='Uniform', FALSE, TRUE), a = a, b = b, alpha = alpha, beta = beta, r = r, lambda = lambda, log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, intercept_move = intercept_move, type_model_prior = type_model_prior, max_parents = max_parents, move_type = move_type, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes)
  # Accept/reject move
  if (min(A,1) >= runif(1)) {
    model <- model_star
    if (is.list(betas_star_included)) {
      betas <- lapply(betas_star_included, function(b_in) replace(numeric(length(model_star)), model_star == 1, b_in))
    } else {
      betas <- replace(numeric(length(model_star)), model_star == 1, betas_star_included)
    }
  }
  
  list(betas = betas, model = model)
}

RJ_MH_Marginal_move = function(Y, X, S = rep(1, length(Y)), model = rep(1,ncol(X)), betas = rep(0,ncol(X)), sigma2 = 1, delta2 = 1, alpha = 0.01, beta = 0.01, a = 0.01, b = 0.01, r = 1, lambda = 1, betas_mean = rep(0,ncol(X)), sigma2_in_betas_prior = FALSE, log = TRUE, intercept_move = FALSE, type_model_prior = 'Uniform', max_parents = ncol(X), allocation_vector = c(1, max(S)), coupling = FALSE){
  # Update Model
  updated_model <- ModelUpdate(k = ncol(X), model = model, max_parents = max_parents, intercept_move = intercept_move)
  model_star <- updated_model$model
  move_type <- updated_model$move
  
  segment_indexes <- if(length(allocation_vector) > 2){
    cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    rep(1, length(Y))
  }
  
  # Construct submatrices for current and proposed models
  X_included <- as.matrix(X[, model == 1])
  X_star <- as.matrix(X[, model_star == 1])
  betas_mean_included <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  betas_mean_star <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model_star == 1])
  A = AcceptanceRatio(Y = Y, X = X_included, X_star = X_star, model = model, model_star = model_star, MCMC_type = 'RJ-Marginal', delta2 = delta2, delta2_star = delta2, betas_0 = betas_mean_included, betas_0_star = betas_mean_star, likelihood = FALSE, marginal = TRUE, model_prior = ifelse(type_model_prior=='Uniform', FALSE, TRUE), beta_prior = FALSE, a = a, b = b, r = r, lambda = lambda, log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, move_type = move_type, intercept_move = intercept_move, type_model_prior = type_model_prior, max_parents = max_parents, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes)
  # Accept/reject move
  if (min(A,1) >= runif(1)) {
    model <- model_star
    
    # Resample betas under accepted model
    betas_star_included <- BetasFCDSample(N = 1, Y = Y, X = X_star, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_mean_star, sigma2_in_betas_prior = TRUE, segment_indexes = segment_indexes)
    
    # Form full betas_star vector (with zeros for excluded covariates)
    if (is.list(betas_star_included)) {
      betas <- lapply(betas_star_included, function(b_in) replace(numeric(length(model_star)), model_star == 1, b_in))
    } else {
      betas <- replace(numeric(length(model_star)), model_star == 1, betas_star_included)
    }
  }
  
  list(betas = betas, model = model)
}

ChangePoint_move <- function(X, Y, S = rep(1, length(Y)), model = rep(1,ncol(X)), betas = rep(0, ncol(X)), sigma2 = 1, delta2 = 1, a = 0.01, b = 0.01, rho = 0.01, iota = 1, betas_mean = rep(0, ncol(X)), type_segment_prior = 'Uniform_card', allocation_vector = c(1, max(S)), max_changepoints = max(S), n_segments = NaN, sigma2_in_betas_prior = TRUE, log = TRUE, coupling = FALSE){
  # Update Changepoints
  params <- UpdateChangepoints(allocation_vector = allocation_vector, max_changepoints = max_changepoints, n_segments = n_segments)
  move_type <- params$move
  allocation_vector_star <- params$allocation_vector
  segment_indexes <- if(length(allocation_vector) > 2){
    cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    rep(1, length(Y))
  }
  #print(segment_indexes)
  segment_indexes_star <- if(length(allocation_vector_star) > 2){
    cut(S, breaks = allocation_vector_star, labels = FALSE, right = FALSE, include.lowest = TRUE)
  } else{
    rep(1, length(Y))
  }
  
  X_included <- X[, model == 1, drop = FALSE]
  betas_included <- lapply(if (is.list(betas)) betas else list(betas), function(b) b[model == 1])
  betas_mean_included <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  betas_mean_star <- lapply(1:max(segment_indexes_star), function(i) {betas_mean_included[[1]]})
  A = AcceptanceRatio(Y = Y, X = X_included, X_star = X_included, model = model, model_star = model, MCMC_type = 'ChangePoint', delta2 = delta2, delta2_star = delta2, betas_0 = betas_mean_included, betas_0_star = betas_mean_star, likelihood = FALSE, marginal = TRUE, beta_prior = FALSE, segment_prior = TRUE, a = a, b = b, rho = rho, iota = iota, log = log, sigma2_in_betas_prior = sigma2_in_betas_prior, move_type = move_type, type_segment_prior = type_segment_prior, max_changepoints = max_changepoints, segment_indexes = segment_indexes, segment_indexes_star = segment_indexes_star, allocation_vector = allocation_vector, allocation_vector_star = allocation_vector_star)
  
  # Accept/reject move
  if (min(A,1) >= runif(1)) {
    allocation_vector <- allocation_vector_star
    # Resample betas under accepted model
    betas_included <- BetasFCDSample(N = 1, Y = Y, X = X_included, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_mean_star, sigma2_in_betas_prior = TRUE, segment_indexes = segment_indexes_star)
    
    # Form full betas_star vector (with zeros for excluded covariates)
    if (is.list(betas_included)) {
      betas <- lapply(betas_included, function(b_in) replace(numeric(length(model)), model == 1, b_in))
    } else {
      betas <- replace(numeric(length(model)), model == 1, betas_included)
    }
    if (is.list(betas_mean)) {
      betas_mean <- lapply(betas_mean_star, function(b_in) replace(numeric(length(model)), model == 1, b_in))
    } else {
      betas_mean <- replace(numeric(length(model)), model == 1, betas_mean_star)
    }
  }
  list(allocation_vector = allocation_vector, betas = betas, betas_mean = betas_mean)
}

RJMCMC = function(data, # data.frame object
                  segmentation = FALSE, # if data segmentation should be applied
                  segmentation_var = 'Time', # variable on which data is segmented
                  n_segments = NaN, # Fixed number of segments (if desired)
                  segmentation_vector = NaN, # Fixed segmentation vector
                  epochs = 1000, # number of epochs of RJMCMC sampling
                  v = 0.2, # parameter for MH moves 
                  p_MCMC = 1, # probability of doing an MCMC move
                  p_RJMCMC = 0.1, # probability of doing an RJMCMC move
                  p_ChangePoint = 0.1, # probability of doing a Change Point Detection move
                  alpha = 0.01, # delta2 scale parameter
                  beta =  0.01, # delta2 rate parameter
                  alpha_xi = 0.01, #xi2 scale parameter
                  beta_xi = 0.01, #xi2 rate parameter
                  a = 0.01, # sigma2 scale parameter
                  b = 0.01, # sigma2 rate parameter
                  r = 1, # number of failures of Neg-Bin model prior
                  lambda = 1, # expected value of number of successes of Poisson/Neg-Bin model prior
                  rho = 0.1, #geometric probability parameter for segment distance prior
                  iota = 1, # expected value of number of successes of Poisson segments' granularity prior
                  intercept_move = FALSE, # if intercept could be deleted too
                  outcome = "Y", # outcome variable
                  covariates = colnames(data)[!colnames(data) %in% outcome], # set of covariates (default: all the variables except outcome)
                  max_parents = length(covariates), # max number of parents 
                  max_changepoints = NaN, #max number of changepoints
                  init_model = character(0), # initial model vector
                  init_betas = rep(0, (length(covariates))), # initial parameters betas
                  betas_mean = rep(0, (length(covariates))), # prior mean of the betas
                  init_sigma2 = 1, # initial sigma2
                  init_delta2 = 1, # initial delta2
                  init_xi2 = 1, # initial xi2
                  init_allocation_vector = NaN, # initial segment set
                  parameters_moves = 'MH', # type of move on the parameter domain ('MH' or 'Gibbs')
                  model_moves = 'MH-Prior', # type of move on the model domain ('MH-Prior', 'MH-FCD' or 'MH-Marginal')
                  type_model_prior = 'Uniform', # type of model prior ('Uniform', 'Uniform_trunc', 'Uniform_card', 'Poisson' or 'Neg-Bin')
                  type_segment_prior = 'Uniform_card', # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
                  type_coupling = 'piece-wise', # type of parameter coupling ('piece-wise' or 'global')
                  log = TRUE, # if logarithm should be applied in acceptance ratio computation
                  coupling = FALSE, # if parameter coupling should be applied
                  sigma2_in_betas_prior = FALSE) { # if sigma2 is present in betas prior
  
  # Experimental Conditions
  X <- as.matrix(data[, (names(data) %in% covariates)])
  Y <- data[, outcome]
  k <- ncol(X)
  
  # Temporal Vector Initialization
  if (segmentation){
    S <- data[, segmentation_var]
  }
  else{
    S <- rep(1, length(Y))
  }
  
  # Changepoint Maximum Number
  max_changepoints <- ifelse(is.nan(max_changepoints), max(unique(data[, segmentation_var])), max_changepoints)
  # Model Initialization
  if (length(init_model) == 0) {
    init_model <- rep(0, k)
    if (max_parents > 0 && max_parents <= k) {
      selected_indices <- sample(k, max_parents)
      init_model[selected_indices] <- 1
    }
  }
  if (intercept_move == FALSE){
    init_model[1] <- 1
    max_parents <- max_parents + 1
  }
  
  # Allocation Vector initialization 
  if (is.nan(n_segments) & any(is.nan(segmentation_vector))){ # unknown changepoints number and changepoints locations
    init_allocation_vector <- c(1,max(S)) # DBN
  }
  else if (!any(is.nan(segmentation_vector))){ # known changepoints number and changepoints locations
    init_allocation_vector <- c(1,segmentation_vector, max(S)) # NH-DBN with fixed changepoints
    n_segments <- length(segmentation_vector) + 1
  }
  else{ # known changepoints number and unknown changepoints locations
    init_allocation_vector <- round(seq(1, max(S), length.out = (n_segments + 1))) # NH-DBN with equi-spaced changepoints
  }
  # Betas Mean Initialization
  betas_mean <- BetasMeanVector(Y = Y, X = X, S = S, model = init_model,
                                 allocation_vector = init_allocation_vector,
                                 init_betas = init_betas,
                                 sigma2 = init_sigma2, delta2 = init_delta2, xi2 = init_xi2, coupled = coupling)

  # Initial Parameters
  model = init_model
  sigma2 = init_sigma2
  betas = list(init_betas) 
  betas = append(betas,replicate(length(init_allocation_vector)-2, rep(0,ncol(X)), simplify = FALSE))
  delta2 = init_delta2
  xi2 = init_xi2
  allocation_vector <- init_allocation_vector
  betas_list = list(betas)
  sigma2s <- list(sigma2)
  delta2s <- list(delta2)
  models <- list(model)
  allocation_vectors <- list(allocation_vector)
  # Sampling Beta parameters
  for (epoch in 1:epochs) {
    # Randomly select between an RJMCMC move and a standard MH move
    if (p_MCMC >= runif(1, min = 0, max = 1)){
      #print('A')
      if (parameters_moves == 'Gibbs'){
        param_set <- Gibbs_move(Y = Y, X = X, S = S, model = model, betas = betas, sigma2 = sigma2, delta2 = delta2, alpha = alpha, beta = beta, a = a, b = b, betas_mean = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, collapsed = ifelse(model_moves == 'Gibbs-Marginal', TRUE, FALSE), allocation_vector = allocation_vector, coupling = coupling)
      }
      else if (parameters_moves == 'MH'){
        param_set <- MH_move(Y = Y, X = X, S = S, model = model, betas = betas, sigma2 = sigma2, delta2 = delta2, alpha = alpha, beta = beta, a = a, b = b, betas_mean = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, log = log, v = v, allocation_vector = allocation_vector, coupling = coupling)
      }
      sigma2 <- param_set$sigma2
      delta2 <- param_set$delta2
      betas <- param_set$betas
    }
    if (p_RJMCMC >= runif(1, min = 0, max = 1)){
      if (model_moves == 'MH-Prior'){
        param_set <- RJ_MH_move(Y = Y, X = X, S = S, model = model, betas = betas, sigma2 = sigma2, delta2 = delta2, alpha = alpha, beta = beta, a = a, b = b, r = r, lambda = lambda, betas_mean = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, log = log, intercept_move = intercept_move, type_model_prior = type_model_prior, max_parents = max_parents, allocation_vector = allocation_vector, coupling = coupling)  
      }
      else if (model_moves == 'MH-FCD'){
        param_set <- RJ_MH_FCD_move(Y = Y, X = X, S = S, model = model, betas = betas, sigma2 = sigma2, delta2 = delta2, alpha = alpha, beta = beta, a = a, b = b, r = r, lambda = lambda, betas_mean = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, log = log, intercept_move = intercept_move, type_model_prior = type_model_prior, max_parents = max_parents, allocation_vector = allocation_vector, coupling = coupling)  
      }
      else if (model_moves == 'MH-Marginal'){
        param_set <- RJ_MH_Marginal_move(Y = Y, X = X, S = S, model = model, betas = betas, sigma2 = sigma2, delta2 = delta2, alpha = alpha, beta = beta, a = a, b = b, r = r, lambda = lambda, betas_mean = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, log = log, intercept_move = intercept_move, type_model_prior = type_model_prior, max_parents = max_parents, allocation_vector = allocation_vector, coupling = coupling)  
      }
      model <- param_set$model
      betas <- param_set$betas
    }
    # if segmentation is not selected or segmentation vector is fixed than the changepoint move is avoided
    if (((segmentation) & (any(is.nan(segmentation_vector)))) & (p_ChangePoint >= runif(1, min = 0, max = 1))){
      #print('C')
      param_set <- ChangePoint_move(Y = Y, X = X, S = S, model = model, betas = betas, sigma2 = sigma2, delta2 = delta2, a = a, b = b, rho = rho, iota = iota, betas_mean = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, log = log, max_changepoints = max_changepoints, allocation_vector = allocation_vector, coupling = coupling, type_segment_prior = type_segment_prior, n_segments = n_segments)  
      allocation_vector <- param_set$allocation_vector
      betas <- param_set$betas
      betas_mean <- param_set$betas_mean
      
    }
    # Updating Markov Chain
    betas_list[[length(betas_list) + 1]] <- betas
    sigma2s[[length(sigma2s) + 1]] <- sigma2
    delta2s[[length(delta2s) + 1]] <- delta2
    models[[length(models) + 1]] <- model
    allocation_vectors[[length(allocation_vectors) + 1]] <- allocation_vector
    #if (epoch %% 1000 == 0) {
    #  cat("Epoch:", epoch, "\n")
    #}
  }
  #cat('\n\n')
  # Return thinned samples post-burn-in
  list(betas = betas_list, sigma2s = sigma2s, delta2s = delta2s, models = models, allocation_vectors = allocation_vectors, covariates = covariates)
}
