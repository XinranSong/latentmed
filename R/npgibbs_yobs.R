## Update for forest_Y, tau2_Y ----

update_forest_Y <- function(Y, XQY, forest_Y) {
  mu_Y <- forest_Y$do_gibbs(XQY, Y, XQY, 1)
  return(list(mu_Y = mu_Y, tau2_Y = forest_Y$get_sigma()^2))
}

update_forest_M <- function(theta_M, XQ, forest_M) {
  mu_theta_M <- forest_M$do_gibbs(XQ, theta_M, XQ, 1)
  return(list(mu_theta_M = mu_theta_M, tau2_M = forest_M$get_sigma()^2))
}

## Update Slopes and Intercepts for M only ----

update_a_b_sigma <- function(theta, Y, sigma2) {
  J <- ncol(Y)
  N <- nrow(Y)
  a <- numeric(J)
  b <- numeric(J)
  a[1] <- 0
  b[1] <- 1
  for (j in 2:J) {
    X_mat <- cbind(1, theta)
    prec_mat <- crossprod(X_mat)/sigma2[j]
    mu_hat <- solve(prec_mat) %*% crossprod(X_mat, Y[, j]/sigma2[j])
    ab_j <- MASS::mvrnorm(1, mu_hat, solve(prec_mat))
    a[j] <- ab_j[1]
    b[j] <- ab_j[2]
  }

  for (j in 1:J) {
    residuals <- Y[, j] - a[j] - b[j] * theta
    sigma2[j] <- 1/rgamma(1, shape = N/2, rate = sum(residuals^2)/2)
  }

  return(list(a = a, b = b, sigma2 = sigma2))
}

update_theta <- function(mu_theta, Y, a, b, sigma2, tau2) {
  N <- nrow(Y)
  J <- ncol(Y)

  Y_centered <- Y - tcrossprod(rep(1, N), a)
  weighted_residuals <- Y_centered * tcrossprod(rep(1, N), b/sigma2)
  data_term <- rowSums(weighted_residuals)

  sum_num <- data_term + mu_theta/tau2
  sum_den <- sum((b^2)/sigma2) + 1/tau2

  post_mean <- sum_num/sum_den
  post_sd <- sqrt(1/sum_den)
  theta_new <- rnorm(N, post_mean, post_sd)

  return(theta_new)
}

update_theta_M <- function(mu_theta_M, theta_M, M, a, b, sigma2_M, tau2_M, Y, forest_Y, XQ) {
  N <- length(theta_M)
  XQY <- cbind(XQ, pnorm(theta_M))
  prior_loglik <- dnorm(x = Y,
                        mean = forest_Y$do_predict(XQY),
                        sd = forest_Y$get_sigma(),
                        log = TRUE)

  theta_M_prop <- update_theta(mu_theta_M, M, a, b, sigma2_M, tau2_M)
  XQY_prop <- cbind(XQ, pnorm(theta_M_prop))
  post_loglik <- dnorm(x = Y,
                       mean = forest_Y$do_predict(XQY_prop),
                       sd = forest_Y$get_sigma(),
                       log = TRUE)
  U <- log(runif(N))
  return(ifelse(U < post_loglik - prior_loglik, theta_M_prop, theta_M))
}

initialize_vars <- function(X, Y, M, n_burn = 2500, n_save = 2500) {

  X1 <- cbind(1, X)
  XQ <- quantile_normalize_bart(X)

  ## Initialize theta_M using factor analysis on M
  fit_theta_M <- gibbs_factor_analysis(Y = M, X = X1, n_iter = n_burn + n_save)
  theta_M <- colMeans(fit_theta_M$theta[-(1:n_burn), ])
  sigma2_M <- colMeans(fit_theta_M$sigma2[-(1:n_burn),])

  ## Initialize regression from X to M
  opts <- Opts(update_s = FALSE)
  hypers_M <- Hypers(X, theta_M)
  forest_M <- MakeForest(hypers = hypers_M, opts = opts)

  mu_theta_M <- forest_M$do_gibbs(XQ, theta_M, XQ, n_burn)
  mu_theta_M <- forest_M$do_gibbs(XQ, theta_M, XQ, n_burn)
  tau2_M <- forest_M$get_sigma()^2

  ## Initialize regression from (X, theta_M) to Y
  hypers_Y <- Hypers(cbind(X, theta_M), Y)
  forest_Y <- MakeForest(hypers = hypers_Y, opts = opts)

  XQY <- cbind(XQ, pnorm(theta_M))
  mu_Y <- forest_Y$do_gibbs(XQY, Y, XQY, n_burn)
  mu_Y <- forest_Y$do_gibbs(XQY, Y, XQY, n_burn)
  tau2_Y <- forest_Y$get_sigma()^2

  ## Initialize slopes and sigma for M
  slopes_M <- update_a_b_sigma(theta_M, M, sigma2_M)
  a_M <- slopes_M$a
  b_M <- slopes_M$b
  sigma2_M <- slopes_M$sigma2

  inits <- list(theta_M = theta_M, forest_Y = forest_Y,
                forest_M = forest_M, a_M = a_M, b_M = b_M,
                sigma2_M = sigma2_M, tau2_Y = tau2_Y,
                tau2_M = tau2_M, XQ = XQ, X = X, Y = Y, M = M)

  return(inits)
}

npgibbs_iter <- function(state) {
  Y <- state$Y
  theta_M <- state$theta_M
  forest_Y <- state$forest_Y
  forest_M <- state$forest_M
  a_M <- state$a_M
  b_M <- state$b_M
  sigma2_M <- state$sigma2_M
  tau2_Y <- state$tau2_Y
  tau2_M <- state$tau2_M
  X <- state$X
  M <- state$M
  XQ <- state$XQ

  XQY <- cbind(XQ, pnorm(theta_M))

  # Update forests
  updated_forest_Y <- update_forest_Y(Y, XQY, forest_Y)
  updated_forest_M <- update_forest_M(theta_M, XQ, forest_M)
  tau2_Y <- updated_forest_Y$tau2_Y
  tau2_M <- updated_forest_M$tau2_M
  mu_theta_M <- updated_forest_M$mu_theta_M

  # Update slopes and sigma for M
  slopes_M <- update_a_b_sigma(theta_M, M, sigma2_M)
  a_M <- slopes_M$a
  b_M <- slopes_M$b
  sigma2_M <- slopes_M$sigma2

  # Update theta_M
  theta_M <- update_theta_M(mu_theta_M = mu_theta_M, theta_M = theta_M, M = M,
                            a = a_M, b = b_M, sigma2_M = sigma2_M,
                            tau2_M = tau2_M, Y = Y,
                            forest_Y = forest_Y, XQ = XQ)

  # Update state
  state$theta_M <- theta_M
  state$forest_Y <- forest_Y
  state$forest_M <- forest_M
  state$a_M <- a_M
  state$b_M <- b_M
  state$sigma2_M <- sigma2_M
  state$tau2_Y <- tau2_Y
  state$tau2_M <- tau2_M

  return(state)
}


## Single-draw Dirichlet sampler for Bayesian bootstrap weights ----
rdirichlet_1 <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

#' @export
run_npgibbs <- function(X, Y, M, treatment_col = 1, n_burn = 2500, n_save = 2500) {
  state <- initialize_vars(X, Y, M, n_burn, n_save)

  n <- nrow(X)
  n_m <- ncol(M)

  samples <- list(
    theta_M = matrix(NA, n_save, n),
    a_M = matrix(NA, n_save, n_m),
    b_M = matrix(NA, n_save, n_m),
    sigma2_M = matrix(NA, n_save, n_m),
    tau2_Y = numeric(n_save),
    tau2_M = numeric(n_save),
    direct_0 = numeric(n_save),
    direct_1 = numeric(n_save),
    indirect_0 = numeric(n_save),
    indirect_1 = numeric(n_save),
    # PATE via Bayesian bootstrap
    direct_0_pate   = numeric(n_save),
    direct_1_pate   = numeric(n_save),
    indirect_0_pate = numeric(n_save),
    indirect_1_pate = numeric(n_save)
  )

  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = n_burn + n_save
  )

  # Burn-in
  for(i in 1:n_burn) {
    pb$tick()
    state <- npgibbs_iter(state)
  }

  # Sampling
  for(i in 1:n_save) {
    pb$tick()
    state <- npgibbs_iter(state)

    samples$theta_M[i,] <- state$theta_M
    samples$a_M[i,] <- state$a_M
    samples$b_M[i,] <- state$b_M
    samples$sigma2_M[i,] <- state$sigma2_M
    samples$tau2_Y[i] <- state$tau2_Y
    samples$tau2_M[i] <- state$tau2_M

    # Counterfactual calculations
    XQ_0 <- XQ_1 <- state$XQ
    XQ_1[,treatment_col] <- 1
    XQ_0[,treatment_col] <- 0

    epsilon_M <- rnorm(n, 0, sqrt(state$tau2_M))
    M_1 <- state$forest_M$do_predict(XQ_1) + epsilon_M
    M_0 <- state$forest_M$do_predict(XQ_0) + epsilon_M

    epsilon_Y <- rnorm(n, 0, sqrt(state$tau2_Y))
    XQY_00 <- cbind(XQ_0, pnorm(M_0))
    XQY_01 <- cbind(XQ_0, pnorm(M_1))
    XQY_10 <- cbind(XQ_1, pnorm(M_0))
    XQY_11 <- cbind(XQ_1, pnorm(M_1))

    Y_00 <- state$forest_Y$do_predict(XQY_00) + epsilon_Y
    Y_01 <- state$forest_Y$do_predict(XQY_01) + epsilon_Y
    Y_10 <- state$forest_Y$do_predict(XQY_10) + epsilon_Y
    Y_11 <- state$forest_Y$do_predict(XQY_11) + epsilon_Y

    samples$direct_0[i] <- mean(Y_10 - Y_00)
    samples$direct_1[i] <- mean(Y_11 - Y_01)
    samples$indirect_0[i] <- mean(Y_01 - Y_00)
    samples$indirect_1[i] <- mean(Y_11 - Y_10)

    ## PATE via Bayesian bootstrap weights
    w <- rdirichlet_1(rep(1, n))

    samples$direct_0_pate[i]   <- sum(w * (Y_10 - Y_00))
    samples$direct_1_pate[i]   <- sum(w * (Y_11 - Y_01))
    samples$indirect_0_pate[i] <- sum(w * (Y_01 - Y_00))
    samples$indirect_1_pate[i] <- sum(w * (Y_11 - Y_10))
  }

  return(samples)
}
