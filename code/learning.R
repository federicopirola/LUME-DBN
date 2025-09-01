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