em_estimate <- function(Y, max_iter=1000, tol=1e-6) {
  # 1. Ensure Y is a matrix (n x J). If Y is a vector, make it a (n x 1) matrix.
  if (is.null(dim(Y))) { 
    Y <- matrix(Y, ncol = 1)
  }
  
  n <- nrow(Y)  # Number of observations
  J <- ncol(Y)  # Number of measures
  
  # 2. Initialize parameters
  mu_theta <- mean(Y[, 1])
  sigma_theta <- var(Y[, 1])
  beta0 <- colMeans(Y)
  beta1 <- rep(1, J)
  sigma_e <- apply(Y, 2, var)
  
  # Pin measure #1
  beta0[1] <- 0
  beta1[1] <- 1
  
  old_ll <- -Inf
  
  for (iter in 1:max_iter) {
    # E-step
    prec <- 1 / sigma_theta + 1 / sigma_e[1]
    if (J > 1) {
      for (j in 2:J) {
        prec <- prec + (beta1[j]^2) / sigma_e[j]
      }
    }
    post_var_theta <- 1 / prec
    
    mean_comp <- mu_theta / sigma_theta + Y[, 1] / sigma_e[1]
    if (J > 1) {
      for (j in 2:J) {
        mean_comp <- mean_comp + beta1[j] * (Y[, j] - beta0[j]) / sigma_e[j]
      }
    }
    post_mean_theta <- post_var_theta * mean_comp
    
    # M-step
    mu_theta <- mean(post_mean_theta)
    sigma_theta <- mean(post_var_theta + (post_mean_theta - mu_theta)^2)
    
    if (J > 1) {
      E_theta <- post_mean_theta
      E_theta2 <- post_var_theta + post_mean_theta^2  
      
      for (j in 2:J) {
        sum_E_theta2 <- sum(E_theta2)
        sum_E_theta <- sum(E_theta)
        sum_Y <- sum(Y[, j])
        sum_YE_theta <- sum(Y[, j] * E_theta)
        
        denom <- n * sum_E_theta2 - (sum_E_theta)^2
        beta1[j] <- (n * sum_YE_theta - sum_E_theta * sum_Y) / denom
        beta0[j] <- (sum_E_theta2 * sum_Y - sum_E_theta * sum_YE_theta) / denom
        
        pred <- beta0[j] + beta1[j] * post_mean_theta
        sigma_e[j] <- mean((Y[, j] - pred)^2 + beta1[j]^2 * post_var_theta)
      }
    }
    
    sigma_e[1] <- mean((Y[, 1] - post_mean_theta)^2 + post_var_theta)
    
    # Compute log-likelihood
    ll <- -n * J / 2 * log(2 * pi)
    ll <- ll - n / 2 * log(sigma_e[1]) - sum((Y[, 1] - post_mean_theta)^2) / (2 * sigma_e[1])
    
    if (J > 1) {
      for (j in 2:J) {
        ll <- ll - n / 2 * log(sigma_e[j])
        pred <- beta0[j] + beta1[j] * post_mean_theta
        ll <- ll - sum((Y[, j] - pred)^2) / (2 * sigma_e[j])
      }
    }
    
    ll <- ll - n / 2 * log(sigma_theta) - sum((post_mean_theta - mu_theta)^2) / (2 * sigma_theta)
    
    if (abs(ll - old_ll) < tol) {
      break
    }
    old_ll <- ll
  }
  
  return(list(
    mu_theta = mu_theta,
    sigma_theta = sigma_theta,
    beta0 = beta0,
    beta1 = beta1,
    sigma_e = sigma_e,
    theta = post_mean_theta,
    iterations = iter,
    loglik = ll
  ))
}



# em_estimate <- function(Y, max_iter=1000, tol=1e-6) {
#   n <- nrow(Y)
#   J <- ncol(Y)
#   
#   # Initialize parameters
#   mu_theta <- mean(Y[,1])
#   sigma_theta <- var(Y[,1])
#   beta0 <- colMeans(Y)
#   beta1 <- rep(1, J)
#   sigma_e <- apply(Y, 2, var)
#   
#   # Fix parameters for first measure
#   beta0[1] <- 0
#   beta1[1] <- 1
#   
#   # Store log-likelihood
#   old_ll <- -Inf
#   
#   for(iter in 1:max_iter) {
#     # E-step
#     prec <- 1/sigma_theta + 1/sigma_e[1]
#     for(j in 2:J) {
#       prec <- prec + (beta1[j]^2)/sigma_e[j]
#     }
#     post_var_theta <- 1/prec
#     
#     mean_comp <- mu_theta/sigma_theta + Y[,1]/sigma_e[1]
#     for(j in 2:J) {
#       mean_comp <- mean_comp + beta1[j] * (Y[,j] - beta0[j])/sigma_e[j]
#     }
#     post_mean_theta <- post_var_theta * mean_comp
#     
#     # M-step
#     # Update mu_theta and sigma_theta
#     mu_theta <- mean(post_mean_theta)
#     sigma_theta <- mean(post_var_theta + (post_mean_theta - mu_theta)^2)
#     
#     # Update parameters for j > 1
#     for(j in 2:J) {
#       # Update beta considering uncertainty in theta
#       E_theta2 <- post_var_theta + post_mean_theta^2
#       E_theta <- post_mean_theta
#       
#       # Compute sufficient statistics
#       sum_E_theta2 <- sum(E_theta2)
#       sum_E_theta <- sum(E_theta)
#       sum_Y <- sum(Y[,j])
#       sum_YE_theta <- sum(Y[,j] * E_theta)
#       
#       # Update beta
#       denom <- n * sum_E_theta2 - sum_E_theta^2
#       beta1[j] <- (n * sum_YE_theta - sum_E_theta * sum_Y) / denom
#       beta0[j] <- (sum_E_theta2 * sum_Y - sum_E_theta * sum_YE_theta) / denom
#       
#       # Update sigma_e
#       pred <- beta0[j] + beta1[j] * post_mean_theta
#       sigma_e[j] <- mean((Y[,j] - pred)^2 + (beta1[j]^2 * post_var_theta))
#     }
#     
#     # Update sigma_e[1]
#     sigma_e[1] <- mean((Y[,1] - post_mean_theta)^2 + post_var_theta)
#     
#     # Compute log-likelihood
#     ll <- -n * J/2 * log(2*pi)
#     ll <- ll - n/2 * log(sigma_e[1]) - sum((Y[,1] - post_mean_theta)^2)/(2*sigma_e[1])
#     for(j in 2:J) {
#       ll <- ll - n/2 * log(sigma_e[j])
#       pred <- beta0[j] + beta1[j] * post_mean_theta
#       ll <- ll - sum((Y[,j] - pred)^2)/(2*sigma_e[j])
#     }
#     ll <- ll - n/2 * log(sigma_theta) - 
#       sum((post_mean_theta - mu_theta)^2)/(2*sigma_theta)
#     
#     # Check convergence
#     if(abs(ll - old_ll) < tol) break
#     old_ll <- ll
#   }
#   
#   return(list(
#     mu_theta=mu_theta,
#     sigma_theta=sigma_theta,
#     beta0=beta0,
#     beta1=beta1,
#     sigma_e=sigma_e,
#     theta=post_mean_theta,
#     iterations=iter,
#     loglik=ll
#   ))
# }