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