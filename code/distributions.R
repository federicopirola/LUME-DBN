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

BetasFCDSample <- function(N = 1, Y, X, segment_indexes, sigma2, delta2, mu_betas = 0, sigma2_in_betas_prior = FALSE) {
  # Loop through segments based on allocation_vector and S
  betas_list <- lapply(1:max(segment_indexes), function(h) {
    X_current <- X[segment_indexes == h, , drop = FALSE]
    Y_current <- Y[segment_indexes == h]
    # Compute the posterior for each segment
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

CumulativeComplementary <- function(X, rho = 0.01, log = TRUE){
  ifelse(log, (X) * log(1-rho), (1-rho)^(X))
}

LogChoose <- function(n, k) {
  lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
}

Likelihood <- function(Y, X, betas = rep(0, ncol(X)), sigma2 = 1, segment_indexes = rep(1, length(Y)), log = TRUE){
  sum(sapply(1:max(segment_indexes), function(i) {MultivariateNormalDensity(X = Y[segment_indexes==i], mu = X[segment_indexes==i, ,drop = FALSE] %*% betas[[i]], Sigma = sigma2, log = log)}))
}

BetaPrior <- function(betas, betas_0 = rep(0, length(betas[[1]])), segment_indexes = rep(1, length(Y)), sigma2 = 1, delta2 = 1, sigma2_in_betas_prior = FALSE, log = TRUE){
  BetaNoise = ifelse(sigma2_in_betas_prior, sigma2*delta2, delta2)
  sum(sapply(1:max(segment_indexes), function(i) {MultivariateNormalDensity(betas[[i]], mu = betas_0[[i]], Sigma = BetaNoise, log = log)}))
}

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