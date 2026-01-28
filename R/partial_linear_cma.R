partial_linear_cma <- function(X, M, Y, num_iter, num_burn, num_trees, priors) {

  N <- nrow(X)
  P <- ncol(M)
  Q <- ncol(Y)


  ##   Initialize parameters
  theta_M <- rnorm(N, mean = 0, sd = 1)
  theta_Y <- rnorm(N, mean = 0, sd = 1)

  ## Set up hyperparameters and options for SoftBART
  opts <- Opts(num_burn = num_burn, num_save = num_iter - num_burn, update_s = FALSE)

  ## Initialize the forests with correct latent variables

  hypers_M <- Hypers(X = X, Y = theta_M, num_tree = num_trees)
  hypers_Y <- Hypers(X = X, Y = theta_Y, num_tree = num_trees)
  f_M_forest <- MakeForest(hypers_M, opts, warn = FALSE)
  f_Y_forest <- MakeForest(hypers_Y, opts, warn = FALSE)


  ## Extract priors
  mu_gamma       <- priors$mu_gamma
  sigma_gamma2   <- priors$sigma_gamma2
  mu_beta_M0     <- priors$mu_beta_M0
  sigma_beta_M02 <- priors$sigma_beta_M02
  mu_beta_M1     <- priors$mu_beta_M1
  sigma_beta_M12 <- priors$sigma_beta_M12
  mu_beta_Y0     <- priors$mu_beta_Y0
  sigma_beta_Y02 <- priors$sigma_beta_Y02
  mu_beta_Y1     <- priors$mu_beta_Y1
  sigma_beta_Y12 <- priors$sigma_beta_Y12
  a_Mp           <- priors$a_Mp
  b_Mp           <- priors$b_Mp
  a_Yq           <- priors$a_Yq
  b_Yq           <- priors$b_Yq
  a_M            <- priors$a_M
  b_M            <- priors$b_M
  a_Y            <- priors$a_Y
  b_Y            <- priors$b_Y

  ## Initialize regression coefficients
  gamma   <- rnorm(1, mean = mu_gamma, sd = sqrt(sigma_gamma2))

  ## ## Initialize f_Y
  f_Y <- rep(0, N)  # Initialize f_Y to zeros of appropriate length

  ## Variances with Inverse-Gamma priors
  sigma_M2   <- 1 / rgamma(1, shape = a_M, rate = b_M)
  sigma_Y2   <- 1 / rgamma(1, shape = a_Y, rate = b_Y)
  sigma_Mp2  <- 1 / rgamma(P, shape = a_Mp, rate = b_Mp)
  sigma_Yq2  <- 1 / rgamma(Q, shape = a_Yq, rate = b_Yq)

  ## Initialize item parameters
  beta_M0 <- rnorm(P, mean = mu_beta_M0, sd = sqrt(sigma_beta_M02))
  beta_M1 <- rnorm(P, mean = mu_beta_M1, sd = sqrt(sigma_beta_M12))
  beta_Y0 <- rnorm(Q, mean = mu_beta_Y0, sd = sqrt(sigma_beta_Y02))
  beta_Y1 <- rnorm(Q, mean = mu_beta_Y1, sd = sqrt(sigma_beta_Y12))

  ## ## Set identifiability constraints for the first items
  beta_M0[1] <- 0
  beta_M1[1] <- 1
  beta_Y0[1] <- 0
  beta_Y1[1] <- 1

  ## Storage for samples after burn-in
  samples <- list(
    theta_Y = matrix(NA, num_iter - num_burn, N),
    theta_M = matrix(NA, num_iter - num_burn, N),
    gamma     = rep(NA, num_iter - num_burn),
    beta_M0   = matrix(NA, num_iter - num_burn, P),
    beta_M1   = matrix(NA, num_iter - num_burn, P),
    beta_Y0   = matrix(NA, num_iter - num_burn, Q),
    beta_Y1   = matrix(NA, num_iter - num_burn, Q),
    sigma_Mp2 = matrix(NA, num_iter - num_burn, P),
    sigma_Yq2 = matrix(NA, num_iter - num_burn, Q),
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
    M_adj      <- t(  (t(M) - beta_M0) * beta_M1 / sigma_Mp2  )
    numerator  <- mu_theta_M / sigma_M2 +
      rowSums(M_adj) +
      gamma * (theta_Y - f_Y) / sigma_Y2
    denom      <- 1 / sigma_M2 + sum(beta_M1^2 / sigma_Mp2) + gamma^2 / sigma_Y2
    mu_post    <- numerator / denom
    sigma_post <- 1 / sqrt(denom)
    theta_M    <- rnorm(n = length(theta_M), mean = mu_post, sd = sigma_post)

    ## 3. Update f_Y using SoftBART
    ## Treat (theta_Y - gamma * theta_M) as response and X
    adjusted_theta_Y <- theta_Y - gamma * theta_M
    f_Y_forest$set_sigma(sqrt(sigma_Y2))
    f_Y <- as.numeric(f_Y_forest$do_gibbs(X, adjusted_theta_Y, X, 1))


    ## 4. Sample theta_Y for each observation
    mu_theta_Y <- f_Y + gamma * theta_M
    Y_adj      <- t((t(Y) - beta_Y0) * beta_Y1 / sigma_Yq2)
    numerator  <- rowSums(Y_adj) + mu_theta_Y / sigma_Y2
    denom      <- (  sum(beta_Y1^2 / sigma_Yq2) + 1 / sigma_Y2  )
    mu_post    <- numerator / denom
    sigma_post <- 1 / sqrt(denom)
    theta_Y    <- rnorm(n = length(theta_Y), mean = mu_post, sd = sigma_post)

    ## 5. Sample gamma
    prior_variance <- Inf
    prior_mean <- 0
    denom <- sum(theta_M^2) / sigma_Y2 + 1 / prior_variance
    numerator <- (prior_mean / prior_variance) + sum(theta_M * (theta_Y - f_Y)) / sigma_Y2
    gamma <- rnorm(1, numerator / denom, 1 / sqrt(denom))

    ## 6. Sample beta_M0 and beta_M1 for mediator items (p >= 2)
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


    ## 7. Sample beta_Y0 and beta_Y1 for outcome items (q >= 2)
    for (q in 2:Q) {
      ## Design matrix for item q
      X_Yq <- cbind(1, theta_Y)
      Y_Yq <- Y[, q]

      ## Prior means and variances
      mu_beta_Yq <- c(mu_beta_Y0, mu_beta_Y1)
      Sigma_beta_Yq_inv <- diag(c(1 / sigma_beta_Y02, 1 / sigma_beta_Y12))

      ## Posterior calculations
      V_beta_Yq <- solve(t(X_Yq) %*% X_Yq / sigma_Yq2[q] + Sigma_beta_Yq_inv)
      m_beta_Yq <- V_beta_Yq %*% (t(X_Yq) %*% Y_Yq / sigma_Yq2[q] +
                                    Sigma_beta_Yq_inv %*% mu_beta_Yq)

      beta_Y_params <- MASS::mvrnorm(1, mu = m_beta_Yq, Sigma = V_beta_Yq)
      beta_Y0[q] <- beta_Y_params[1]
      beta_Y1[q] <- beta_Y_params[2]
    }

    ## 8. Sample sigma_Mp2 (Variances for Mediator Items)
    for (p in 1:P) {
      shape_Mp <- a_Mp + N / 2
      rate_Mp <- b_Mp + 0.5 * sum((M[, p] - beta_M0[p] - beta_M1[p] * theta_M)^2)
      sigma_Mp2[p] <- 1 / rgamma(1, shape = shape_Mp, rate = rate_Mp)
    }


    ## 9. Sample sigma_Yq2 (Variances for Outcome Items)
    for (q in 1:Q) {
      shape_Yq <- a_Yq + N / 2
      rate_Yq <- b_Yq + 0.5 * sum((Y[, q] - beta_Y0[q] - beta_Y1[q] * theta_Y)^2)
      sigma_Yq2[q] <- 1 / rgamma(1, shape = shape_Yq, rate = rate_Yq)
    }

    ## 10. Sample sigma_M2 (Variance for Mediator Latent Variable)
    residuals_M <- theta_M - f_M
    shape_M <- a_M + N / 2
    rate_M <- b_M + 0.5 * sum(residuals_M^2)
    sigma_M2 <- 1 / rgamma(1, shape = shape_M, rate = rate_M)

    ## 11. Sample sigma_Y2 (Variance for Outcome Latent Variable)
    residuals_Y <- theta_Y - f_Y - gamma * theta_M
    shape_Y <- a_Y + N / 2
    rate_Y <- b_Y + 0.5 * sum(residuals_Y^2)
    sigma_Y2 <- 1 / rgamma(1, shape = shape_Y, rate = rate_Y)

    ## Direct and Indirect Effects Calculation
    ## Compute f_M(Z = 1, W) and f_M(Z = 0, W)
    f_M_1 <- f_M_forest$do_predict(cbind(1, X[,-1]))  # Z = 1
    f_M_0 <- f_M_forest$do_predict(cbind(0, X[,-1]))  #  Z = 0

    ## Compute f_Y(Z = 1, W) and f_Y(Z = 0, W)
    f_Y_1 <- f_Y_forest$do_predict(cbind(1, X[,-1]))  #  Z = 1
    f_Y_0 <- f_Y_forest$do_predict(cbind(0, X[,-1]))  # Z = 0

    zeta <- mean(f_Y_1 - f_Y_0) ## direct
    delta <- mean(gamma * (f_M_1 - f_M_0)) ##indirect
    overall <- zeta + delta

    ## Store samples after burn-in period
    if (iter > num_burn) {
      index <- iter - num_burn
      samples$theta_M[index, ]   <- theta_M
      samples$theta_Y[index, ]   <- theta_Y

      samples$gamma[index]       <- gamma
      samples$beta_M0[index, ]   <- beta_M0
      samples$beta_M1[index, ]   <- beta_M1
      samples$beta_Y0[index, ]   <- beta_Y0
      samples$beta_Y1[index, ]   <- beta_Y1
      samples$sigma_Mp2[index, ] <- sigma_Mp2
      samples$sigma_Yq2[index, ] <- sigma_Yq2
      samples$sigma_M2[index]    <- sigma_M2
      samples$sigma_Y2[index]    <- sigma_Y2
      samples$zeta[index] <- zeta
      samples$delta[index] <- delta
      samples$overall[index] <- overall
    }
  }

  return(samples)

}
