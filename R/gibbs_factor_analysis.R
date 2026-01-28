gibbs_factor_analysis <- function(Y, X, n_iter = 1000) {
  N <- nrow(Y)
  J <- ncol(Y)
  P <- ncol(X)
  
  # Initialize parameters with corrected theta initialization
  theta <- Y[, 1] + rnorm(N, 0, sd = 0.1 * sd(Y[, 1]))
  
  # Initialize regression parameters
  gamma <- coef(lm(theta ~ X - 1))
  tau2 <- var(theta - as.numeric(X %*% gamma))
  
  a <- numeric(J)
  a[1] <- 0
  b <- numeric(J)
  b[1] <- 1
  sigma2 <- numeric(J)
  sigma2[1] <- var(Y[, 1] - theta)
  
  # Initialize other parameters
  for (j in 2:J) {
    lm_j <- lm(Y[, j] ~ theta)
    a[j] <- coef(lm_j)[1]
    b[j] <- coef(lm_j)[2]
    sigma2[j] <- var(resid(lm_j))
  }
  
  # Storage for samples
  samples <- list(
    theta = matrix(NA, nrow = n_iter, ncol = N),
    a = matrix(NA, nrow = n_iter, ncol = J),
    b = matrix(NA, nrow = n_iter, ncol = J),
    sigma2 = matrix(NA, nrow = n_iter, ncol = J),
    gamma = matrix(NA, nrow = n_iter, ncol = P),
    tau2 = numeric(n_iter)
  )
  
  # Gibbs sampling loop
  for (iter in 1:n_iter) {
    # 1. Sample theta_i
    for (i in 1:N) {
      prior_mean <- X[i, ] %*% gamma
      sum_num <- sum((Y[i, ] - a) * b / sigma2) + prior_mean/tau2
      sum_den <- sum(b^2 / sigma2) + 1/tau2
      post_var <- 1/sum_den
      post_mean <- sum_num * post_var
      theta[i] <- rnorm(1, post_mean, sqrt(post_var))
    }
    
    # 2. Sample a_j and b_j (j > 1)
    for (j in 2:J) {
      X_mat <- cbind(1, theta)
      prec_mat <- crossprod(X_mat)/sigma2[j]
      mu_hat <- solve(prec_mat) %*% crossprod(X_mat, Y[, j]/sigma2[j])
      ab_j <- MASS::mvrnorm(1, mu_hat, solve(prec_mat))
      a[j] <- ab_j[1]
      b[j] <- ab_j[2]
    }
    
    # 3. Sample sigma2_j
    for (j in 1:J) {
      residuals <- Y[, j] - a[j] - b[j] * theta
      sigma2[j] <- 1/rgamma(1, shape = N/2, rate = sum(residuals^2)/2)
    }
    
    # 4. Sample gamma
    XtX <- crossprod(X)
    Sigma_gamma <- tau2 * solve(XtX + diag(1e-6, P))  # Regularization for stability
    mu_gamma <- solve(XtX + diag(1e-6, P)) %*% crossprod(X, theta)
    gamma <- MASS::mvrnorm(1, mu_gamma, Sigma_gamma)
    
    # 5. Sample tau2
    residuals <- theta - X %*% gamma
    tau2 <- 1/rgamma(1, shape = N/2, rate = sum(residuals^2)/2)
    
    # Store samples
    samples$theta[iter, ] <- theta
    samples$a[iter, ] <- a
    samples$b[iter, ] <- b
    samples$sigma2[iter, ] <- sigma2
    samples$gamma[iter, ] <- gamma
    samples$tau2[iter] <- tau2
    
    if(iter %% 100 == 0) {
      cat("Finishint iteration ", iter, "\n")
    }
    
  }
  
  return(samples)
}

# Validation script
# set.seed(12)
# N <- 1000
# J <- 3
# 
# # Generate covariates (including intercept)
# X <- cbind(1, rnorm(N), rnorm(N))  # 2 covariates + intercept
# true_gamma <- c(0.5, -0.3, 0.2)    # True regression coefficients
# true_tau2 <- 0.8
# true_a <- c(0, 1, -1)
# true_b <- c(1, 0.5, 1.5)
# true_sigma2 <- c(0.5, 1, 0.8)
# 
# # Generate data
# theta_true <- X %*% true_gamma + rnorm(N, 0, sqrt(true_tau2))
# Y <- matrix(0, N, J)
# for (j in 1:J) {
#   Y[, j] <- true_a[j] + true_b[j] * theta_true + rnorm(N, 0, sqrt(true_sigma2[j]))
# }
# 
# # Run sampler
# samples <- gibbs_factor_analysis(Y, X, n_iter = 3000)
# 
# # Analysis of results
# burnin <- 500
# post_gamma <- colMeans(samples$gamma[-(1:burnin),])
# post_tau2 <- mean(samples$tau2[-(1:burnin)])
# 
# cat("True gamma:", true_gamma, "\n")
# cat("Estimated gamma:", post_gamma, "\n\n")
# cat("True tau2:", true_tau2, "\n")
# cat("Estimated tau2:", post_tau2, "\n")