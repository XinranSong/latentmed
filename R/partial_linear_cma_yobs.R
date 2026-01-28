partial_linear_cma_yobs <- function(X, M, Y, num_iter, num_burn, num_trees, priors) {
  
  N <- nrow(X)
  P <- ncol(M)
  
  ## Initialize parameters
  theta_M <- rnorm(N, mean = 0, sd = 1)
  
  ## Set up hyperparameters and options for SoftBART
  opts <- Opts(num_burn = num_burn, num_save = num_iter - num_burn, update_s = FALSE)  
  
  ## Initialize the forests with correct latent variables
  hypers_M <- Hypers(X = X, Y = theta_M, num_tree = num_trees)
  hypers_Y <- Hypers(X = X, Y = Y, num_tree = num_trees)
  f_M_forest <- MakeForest(hypers_M, opts, warn = FALSE)
  f_Y_forest <- MakeForest(hypers_Y, opts, warn = FALSE)
  
  ## Extract priors
  mu_gamma       <- priors$mu_gamma
  sigma_gamma2   <- priors$sigma_gamma2
  mu_beta_M0     <- priors$mu_beta_M0
  sigma_beta_M02 <- priors$sigma_beta_M02
  mu_beta_M1     <- priors$mu_beta_M1
  sigma_beta_M12 <- priors$sigma_beta_M12
  a_Mp           <- priors$a_Mp
  b_Mp           <- priors$b_Mp
  a_M            <- priors$a_M
  b_M            <- priors$b_M
  a_Y            <- priors$a_Y
  b_Y            <- priors$b_Y
  
  ## Initialize regression coefficients
  gamma <- rnorm(1, mean = mu_gamma, sd = sqrt(sigma_gamma2))
  
  ## Initialize f_Y and f_M
  f_Y <- rep(0, N)
  f_M <- rep(0, N)
  
  ## Variances with Inverse-Gamma priors
  sigma_M2   <- 1 / rgamma(1, shape = a_M, rate = b_M)
  sigma_Y2   <- 1 / rgamma(1, shape = a_Y, rate = b_Y)
  sigma_Mp2  <- 1 / rgamma(P, shape = a_Mp, rate = b_Mp)
  
  ## Initialize item parameters for M
  beta_M0 <- rnorm(P, mean = mu_beta_M0, sd = sqrt(sigma_beta_M02))
  beta_M1 <- rnorm(P, mean = mu_beta_M1, sd = sqrt(sigma_beta_M12))
  
  ## Set identifiability constraints for the first item of M
  beta_M0[1] <- 0
  beta_M1[1] <- 1
  
  ## Storage for samples after burn-in
  samples <- list(
    theta_M = matrix(NA, num_iter - num_burn, N),
    gamma     = rep(NA, num_iter - num_burn),
    beta_M0   = matrix(NA, num_iter - num_burn, P),
    beta_M1   = matrix(NA, num_iter - num_burn, P),
    sigma_Mp2 = matrix(NA, num_iter - num_burn, P),
    sigma_M2  = rep(NA, num_iter - num_burn),
    sigma_Y2  = rep(NA, num_iter - num_burn),
    zeta = rep(NA, num_iter - num_burn),
    delta = rep(NA, num_iter - num_burn),
    overall = rep(NA, num_iter - num_burn)
  )
  
  ## Gibbs Sampling Loop
  pb <- progress_bar$new(
    format = "  Sampling [:bar] :percent eta: :eta",
    total = num_iter, clear = FALSE, width= 60)
  
  for (iter in 1:num_iter) {
    pb$tick()
    
    ## 1. Update f_M using SoftBART
    ## Treat theta_M as response and X as predictors
    f_M_forest$set_sigma(sqrt(sigma_M2))
    f_M <- f_M_forest$do_gibbs(X, theta_M, X, 1)  
    
    ## 2. Sample theta_M for each observation
    mu_theta_M <- f_M
    M_adj      <- t((t(M) - beta_M0) * beta_M1 / sigma_Mp2)
    numerator  <- mu_theta_M / sigma_M2 +
      rowSums(M_adj) +
      gamma * (Y - f_Y) / sigma_Y2
    denom      <- 1 / sigma_M2 + sum(beta_M1^2 / sigma_Mp2) + gamma^2 / sigma_Y2
    mu_post    <- numerator / denom
    sigma_post <- 1 / sqrt(denom)
    theta_M    <- rnorm(n = length(theta_M), mean = mu_post, sd = sigma_post)
    
    ## 3. Update f_Y using SoftBART
    ## Treat (Y - gamma * theta_M) as response
    adjusted_Y <- Y - gamma * theta_M
    f_Y_forest$set_sigma(sqrt(sigma_Y2))
    f_Y <- as.numeric(f_Y_forest$do_gibbs(X, adjusted_Y, X, 1))
    
    ## 4. Sample gamma - relationship between theta_M and Y
    prior_variance <- sigma_gamma2
    prior_mean <- mu_gamma
    denom <- sum(theta_M^2) / sigma_Y2 + 1 / prior_variance
    numerator <- (prior_mean / prior_variance) + sum(theta_M * (Y - f_Y)) / sigma_Y2
    gamma <- rnorm(1, numerator / denom, 1 / sqrt(denom))
    
    ## 5. Sample beta_M0 and beta_M1 for mediator items (p >= 2)
    for (p in 2:P) {
      ## Design matrix for item p
      X_Mp <- cbind(1, theta_M)
      Y_Mp <- M[, p]
      
      ## Prior means and variances
      mu_beta_Mp <- c(mu_beta_M0, mu_beta_M1)
      Sigma_beta_Mp_inv <- diag(c(1 / sigma_beta_M02, 1 / sigma_beta_M12))
      
      ## Posterior calculations
      V_beta_Mp <- solve(t(X_Mp) %*% X_Mp / sigma_Mp2[p] + Sigma_beta_Mp_inv)
      m_beta_Mp <- V_beta_Mp %*% (t(X_Mp) %*% Y_Mp / sigma_Mp2[p] +
                                    Sigma_beta_Mp_inv %*% mu_beta_Mp)
      
      beta_M_params <- MASS::mvrnorm(1, mu = m_beta_Mp, Sigma = V_beta_Mp)
      beta_M0[p] <- beta_M_params[1]
      beta_M1[p] <- beta_M_params[2]
    }
    
    ## 6. Sample sigma_Mp2 (Variances for Mediator Items)
    for (p in 1:P) {
      shape_Mp <- a_Mp + N / 2
      rate_Mp <- b_Mp + 0.5 * sum((M[, p] - beta_M0[p] - beta_M1[p] * theta_M)^2)
      sigma_Mp2[p] <- 1 / rgamma(1, shape = shape_Mp, rate = rate_Mp)
    }
    
    ## 7. Sample sigma_M2 (Variance for Mediator Latent Variable)
    residuals_M <- theta_M - f_M
    shape_M <- a_M + N / 2
    rate_M <- b_M + 0.5 * sum(residuals_M^2)
    sigma_M2 <- 1 / rgamma(1, shape = shape_M, rate = rate_M)
    
    ## 8. Sample sigma_Y2 (Variance for Outcome)
    residuals_Y <- Y - f_Y - gamma * theta_M
    shape_Y <- a_Y + N / 2
    rate_Y <- b_Y + 0.5 * sum(residuals_Y^2)
    sigma_Y2 <- 1 / rgamma(1, shape = shape_Y, rate = rate_Y)
    
    ## Direct and Indirect Effects Calculation
    ## Compute f_M(Z = 1, W) and f_M(Z = 0, W)
    f_M_1 <- f_M_forest$do_predict(cbind(1, X[,-1]))  # Z = 1
    f_M_0 <- f_M_forest$do_predict(cbind(0, X[,-1]))  # Z = 0
    
    ## Compute f_Y(Z = 1, W) and f_Y(Z = 0, W)
    f_Y_1 <- f_Y_forest$do_predict(cbind(1, X[,-1]))  # Z = 1
    f_Y_0 <- f_Y_forest$do_predict(cbind(0, X[,-1]))  # Z = 0
    
    zeta <- mean(f_Y_1 - f_Y_0)      ## direct effect
    delta <- mean(gamma * (f_M_1 - f_M_0))  ## indirect effect
    overall <- zeta + delta          ## total effect
    
    ## Store samples after burn-in period
    if (iter > num_burn) {
      index <- iter - num_burn
      samples$theta_M[index, ] <- theta_M
      samples$gamma[index]     <- gamma
      samples$beta_M0[index, ] <- beta_M0
      samples$beta_M1[index, ] <- beta_M1
      samples$sigma_Mp2[index, ] <- sigma_Mp2
      samples$sigma_M2[index]  <- sigma_M2
      samples$sigma_Y2[index]  <- sigma_Y2
      samples$zeta[index]      <- zeta
      samples$delta[index]     <- delta
      samples$overall[index]   <- overall
    }
  }
  
  return(samples)
}