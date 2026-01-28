## Update for forest_Y, tau2_Y ----

update_forest_Y <- function(Y, XQY, forest_Y) {
  mu_Y <- forest_Y$do_gibbs(XQY, Y, XQY, 1)
  return(list(mu_Y = mu_Y, tau2_Y = forest_Y$get_sigma()^2))
}

update_forest_M <- function(theta_M, XQ, forest_M) {
  mu_theta_M <- forest_M$do_gibbs(XQ, theta_M, XQ, 1)
  return(list(mu_theta_M = mu_theta_M, tau2_M = forest_M$get_sigma()^2))
}

## Updated Slopes and Intercepts (Generalized for M or Y) ----

update_a_b_sigma <- function(theta, Data_Mat, sigma2) {
  Data_Mat <- as.matrix(Data_Mat)
  J <- ncol(Data_Mat)
  N <- nrow(Data_Mat)
  a <- numeric(J)
  b <- numeric(J)
  a[1] <- 0
  b[1] <- 1

  for (j in 2:J) {
    X_mat <- cbind(1, theta)
    prec_mat <- (crossprod(X_mat) + diag(1e-6, 2)) / sigma2[j]
    rhs <- crossprod(X_mat, Data_Mat[, j]) / sigma2[j]
    mu_hat <- solve(prec_mat, rhs)
    ab_j <- MASS::mvrnorm(1, mu_hat, solve(prec_mat))
    a[j] <- ab_j[1]
    b[j] <- ab_j[2]
  }

  for (j in 1:J) {
    residuals <- Data_Mat[, j] - a[j] - b[j] * theta
    sigma2[j] <- 1 / rgamma(1, shape = N/2, rate = sum(residuals^2)/2)
  }

  return(list(a = a, b = b, sigma2 = sigma2))
}

## Updated Latent Score Update ----

update_theta_latent <- function(mu_theta, Data_Mat, a, b, sigma2, tau2) {
  Data_Mat <- as.matrix(Data_Mat)
  N <- nrow(Data_Mat)
  Data_centered <- Data_Mat - tcrossprod(rep(1, N), a)
  weighted_residuals <- Data_centered * tcrossprod(rep(1, N), b/sigma2)
  data_term <- rowSums(weighted_residuals)

  sum_num <- data_term + mu_theta/tau2
  sum_den <- sum((b^2)/sigma2) + 1/tau2

  post_mean <- sum_num/sum_den
  post_sd <- sqrt(1/sum_den)
  return(rnorm(N, post_mean, post_sd))
}

## Updated Metropolis-Hastings for Theta_M ----

update_theta_M <- function(mu_theta_M, theta_M, M, a_M, b_M, sigma2_M, tau2_M,
                           Y, a_Y, b_Y, sigma2_Y, forest_Y, XQ, Y_is_latent) {
  N <- length(theta_M)

  # Current Likelihood
  XQY_curr <- cbind(XQ, pnorm(theta_M))
  pred_curr <- forest_Y$do_predict(XQY_curr)

  # loglik under current theta_M
  if (!Y_is_latent) {
    # Y is observed vector; forest_Y models Y directly
    ll_cur <- dnorm(Y, mean = pred_curr, sd = forest_Y$get_sigma(), log = TRUE)
  } else {
    # Y is indicators matrix; forest_Y models theta_Y, and Y|theta_Y uses measurement model
    # pred_curr is theta_Y mean
    mu_mat <- tcrossprod(pred_curr, b_Y) + tcrossprod(rep(1, N), a_Y)  # N x Jy
    sd_mat <- matrix(sqrt(sigma2_Y), nrow = N, ncol = length(sigma2_Y), byrow = TRUE)
    ll_cur <- rowSums(dnorm(Y, mean = mu_mat, sd = sd_mat, log = TRUE))
  }

  # propose theta_M from measurement-model conditional
  theta_prop <- update_theta_latent(mu_theta_M, M, a_M, b_M, sigma2_M, tau2_M)

  XQY_prop <- cbind(XQ, pnorm(theta_prop))
  pred_prop <- forest_Y$do_predict(XQY_prop)

  if (!Y_is_latent) {
    ll_prop <- dnorm(Y, mean = pred_prop, sd = forest_Y$get_sigma(), log = TRUE)
  } else {
    mu_mat_p <- tcrossprod(pred_prop, b_Y) + tcrossprod(rep(1, N), a_Y)
    sd_mat_p <- matrix(sqrt(sigma2_Y), nrow = N, ncol = length(sigma2_Y), byrow = TRUE)
    ll_prop <- rowSums(dnorm(Y, mean = mu_mat_p, sd = sd_mat_p, log = TRUE))
  }

  U <- log(runif(N))
  accept <- (U < (ll_prop - ll_cur))
  return(ifelse(accept, theta_prop, theta_M))
}



initialize_vars <- function(X, Y, M, n_burn = 2500, n_save = 2500, forest_warm = 100) {

  X1 <- cbind(1, X)
  XQ <- quantile_normalize_bart(X)
  Y_mat <- as.matrix(Y); M_mat <- as.matrix(M)

  # 2x2 Identification Check

  M_is_latent <- ncol(M_mat) >= 3
  Y_is_latent <- ncol(Y_mat) >= 3
  if(ncol(M_mat) == 2) stop("M has 2 columns. Need 1 or >=3.")
  if(ncol(Y_mat) == 2) stop("Y has 2 columns. Need 1 or >=3.")


  ## ---- theta_M init
  if (M_is_latent) {
    fitM <- gibbs_factor_analysis(Y = M_mat, X = X1, n_iter = n_burn + n_save)
    theta_M <- colMeans(fitM$theta[(n_burn + 1):(n_burn + n_save), , drop = FALSE])
    sigma2_M <- colMeans(fitM$sigma2[(n_burn + 1):(n_burn + n_save), , drop = FALSE])
  } else {
    theta_M <- as.numeric(M_mat[, 1])
    sigma2_M <- 1
  }

  ## ---- theta_Y init
  if (Y_is_latent) {
    fitY <- gibbs_factor_analysis(Y = Y_mat, X = X1, n_iter = n_burn + n_save)
    theta_Y <- colMeans(fitY$theta[(n_burn + 1):(n_burn + n_save), , drop = FALSE])
    sigma2_Y <- colMeans(fitY$sigma2[(n_burn + 1):(n_burn + n_save), , drop = FALSE])
  } else {
    theta_Y <- as.numeric(Y_mat[, 1])
    sigma2_Y <- 1
  }

  ## ---- forests
  opts <- Opts(update_s = FALSE)

  forest_M <- MakeForest(hypers = Hypers(X, theta_M), opts = opts)

  # forest_Y always models the "target" used in updates:
  # - observed Y => theta_Y is Y vector
  # - latent Y => theta_Y is latent score
  forest_Y <- MakeForest(hypers = Hypers(cbind(X, theta_M), theta_Y), opts = opts)

  ## ---- warm-start forests so sigmas are meaningful
  forest_M$do_gibbs(XQ, theta_M, XQ, forest_warm)

  XQY <- cbind(XQ, pnorm(theta_M))
  forest_Y$do_gibbs(XQY, theta_Y, XQY, forest_warm)

  ## ---- measurement params
  if (M_is_latent) {
    slopes_M <- update_a_b_sigma(theta_M, M_mat, sigma2_M)
    a_M <- slopes_M$a; b_M <- slopes_M$b; sigma2_M <- slopes_M$sigma2
  } else {
    a_M <- 0; b_M <- 1; sigma2_M <- 1
  }

  if (Y_is_latent) {
    slopes_Y <- update_a_b_sigma(theta_Y, Y_mat, sigma2_Y)
    a_Y <- slopes_Y$a; b_Y <- slopes_Y$b; sigma2_Y <- slopes_Y$sigma2
  } else {
    a_Y <- 0; b_Y <- 1; sigma2_Y <- 1
  }

  list(
    X = X, XQ = XQ,
    Y = Y_mat, M = M_mat,
    theta_M = theta_M,
    theta_Y = theta_Y,
    forest_M = forest_M,
    forest_Y = forest_Y,
    a_M = a_M, b_M = b_M, sigma2_M = sigma2_M,
    a_Y = a_Y, b_Y = b_Y, sigma2_Y = sigma2_Y,
    tau2_M = forest_M$get_sigma()^2,
    tau2_Y = forest_Y$get_sigma()^2,
    mu_theta_M = forest_M$do_predict(XQ),
    mu_theta_Y = forest_Y$do_predict(XQY),
    M_is_latent = M_is_latent,
    Y_is_latent = Y_is_latent
  )
}

## ---------------------------
## One iteration
## ---------------------------

npgibbs_iter <- function(state) {

  XQ <- state$XQ


  # update forests
  up_M <- update_forest_M(state$theta_M, XQ, state$forest_M)
  state$mu_theta_M <- up_M$mu_theta_M
  state$tau2_M <- up_M$tau2_M

  XQY <- cbind(XQ, pnorm(state$theta_M))
  up_Y <- update_forest_Y(state$theta_Y, XQY, state$forest_Y)
  state$mu_theta_Y <- up_Y$mu_Y
  state$tau2_Y <- up_Y$tau2_Y

  # update measurement model + latent scores
  if (state$M_is_latent) {
    sl_M <- update_a_b_sigma(state$theta_M, state$M, state$sigma2_M)
    state$a_M <- sl_M$a; state$b_M <- sl_M$b; state$sigma2_M <- sl_M$sigma2
    # theta_M updated with MH step that uses Y side likelihood
    state$theta_M <- update_theta_M(
      mu_theta_M = state$mu_theta_M,
      theta_M    = state$theta_M,
      M      = state$M,
      a_M        = state$a_M,
      b_M        = state$b_M,
      sigma2_M   = state$sigma2_M,
      tau2_M     = state$tau2_M,
      Y      = if (state$Y_is_latent) state$Y else as.numeric(state$Y[, 1]),
      a_Y        = state$a_Y,
      b_Y        = state$b_Y,
      sigma2_Y   = state$sigma2_Y,
      forest_Y   = state$forest_Y,
      XQ         = XQ,
      Y_is_latent = state$Y_is_latent
    )
  }

  if (state$Y_is_latent) {
    sl_Y <- update_a_b_sigma(state$theta_Y, state$Y, state$sigma2_Y)
    state$a_Y <- sl_Y$a; state$b_Y <- sl_Y$b; state$sigma2_Y <- sl_Y$sigma2
    # update theta_Y via conjugate (given mu_theta_Y from forest_Y)
    state$theta_Y <- update_theta_latent(state$mu_theta_Y, state$Y,
                                         state$a_Y, state$b_Y, state$sigma2_Y,
                                         state$tau2_Y)
  }

  state
}

## ---------------------------
## Dirichlet weights (Bayesian bootstrap)
## ---------------------------
rdirichlet_1 <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

## ---------------------------
## Main runner
## ---------------------------

#' @export
run_npgibbs_general <- function(X, Y, M,
                                treatment_col = 1,
                                n_burn = 2500,
                                n_save = 2500,
                                forest_warm = 100) {

  state <- initialize_vars(X, Y, M, n_burn = n_burn, n_save = n_save, forest_warm = forest_warm)

  n <- nrow(state$X)

  # outcome used for effects:
  # - if Y observed: effects on observed Y
  # - if Y latent: effects on theta_Y
  y_is_latent <- state$Y_is_latent

  # allocate samples
  samples <- list(
    theta_M = matrix(NA, n_save, n),
    theta_Y = if (y_is_latent) matrix(NA, n_save, n) else NULL,
    tau2_Y = numeric(n_save),
    tau2_M = numeric(n_save),
    direct_0 = numeric(n_save),
    direct_1 = numeric(n_save),
    indirect_0 = numeric(n_save),
    indirect_1 = numeric(n_save),
    direct_0_pate = numeric(n_save),
    direct_1_pate = numeric(n_save),
    indirect_0_pate = numeric(n_save),
    indirect_1_pate = numeric(n_save)
  )

  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = n_burn + n_save
  )

  # burn-in
  for (i in 1:n_burn) {
    state <- npgibbs_iter(state)
    pb$tick()
  }

  # sampling
  for (i in 1:n_save) {
    state <- npgibbs_iter(state)
    pb$tick()

    samples$theta_M[i, ] <- state$theta_M
    if (y_is_latent) samples$theta_Y[i, ] <- state$theta_Y
    samples$tau2_Y[i] <- state$tau2_Y
    samples$tau2_M[i] <- state$tau2_M

    # Counterfactual calculations
    XQ_0 <- XQ_1 <- state$XQ
    XQ_1[, treatment_col] <- 1
    XQ_0[, treatment_col] <- 0

    eps_M <- rnorm(n, 0, sqrt(state$tau2_M))
    M_1 <- state$forest_M$do_predict(XQ_1) + eps_M
    M_0 <- state$forest_M$do_predict(XQ_0) + eps_M

    eps_Y <- rnorm(n, 0, sqrt(state$tau2_Y))

    XQY_00 <- cbind(XQ_0, pnorm(M_0))
    XQY_01 <- cbind(XQ_0, pnorm(M_1))
    XQY_10 <- cbind(XQ_1, pnorm(M_0))
    XQY_11 <- cbind(XQ_1, pnorm(M_1))

    # predicted outcome target (observed Y or theta_Y)
    Yt_00 <- state$forest_Y$do_predict(XQY_00) + eps_Y
    Yt_01 <- state$forest_Y$do_predict(XQY_01) + eps_Y
    Yt_10 <- state$forest_Y$do_predict(XQY_10) + eps_Y
    Yt_11 <- state$forest_Y$do_predict(XQY_11) + eps_Y

    samples$direct_0[i]   <- mean(Yt_10 - Yt_00)
    samples$direct_1[i]   <- mean(Yt_11 - Yt_01)
    samples$indirect_0[i] <- mean(Yt_01 - Yt_00)
    samples$indirect_1[i] <- mean(Yt_11 - Yt_10)

    # PATE via Bayesian bootstrap
    w <- rdirichlet_1(rep(1, n))
    samples$direct_0_pate[i]   <- sum(w * (Yt_10 - Yt_00))
    samples$direct_1_pate[i]   <- sum(w * (Yt_11 - Yt_01))
    samples$indirect_0_pate[i] <- sum(w * (Yt_01 - Yt_00))
    samples$indirect_1_pate[i] <- sum(w * (Yt_11 - Yt_10))
  }

  samples
}
