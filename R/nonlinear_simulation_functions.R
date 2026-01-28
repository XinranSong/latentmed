#' @export
generate_data_nonlinear <- function(sigma_factor = 1, nonlinear_strength = 5, N_factor = 1) {
  # Step 1: Estimate latent variables using factor analysis
  fa_result <- factanal(M, factors = 1)
  loadings <- fa_result$loadings[,1]
  factor_scores <- factanal(M, factors = 1, scores = "regression")$scores[,1]

  # Step 2: Fit models to estimate parameters
  composite_model <- data.frame(Z = X[,2], X.1 = X[,3],
                                X.2 = X[,4], X.3 = X[,5],
                                f = factor_scores)
  theta_model <- lm(f ~ (.)^2, data = composite_model)

  # Create orthogonal nonlinear effect using X.2
  X2_squared <- X[,4]^2

  # Residualize X.2 squared against X.2 to ensure orthogonality
  X2_squared_model <- lm(X2_squared ~ X[,4])
  X2_squared_orthogonal <- residuals(X2_squared_model)

  # Standardize the orthogonal nonlinear effect
  X2_squared_orthogonal_std <- scale(X2_squared_orthogonal)[,1]

  # Item models remain unchanged
  item_params <- list()
  for(p in 1:ncol(M)) {
    model <- lm(M[,p] ~ factor_scores)
    if(p == 1) {
      item_params[[p]] <- list(
        intercept = 0,
        loading = 1,
        residual_sd = sd(residuals(model))
      )
    } else {
      item_params[[p]] <- list(
        intercept = coef(model)[1],
        loading = coef(model)[2],
        residual_sd = sd(residuals(model))
      )
    }
  }

  # Modify outcome model to include the standardized orthogonal nonlinear effect
  outcome_df <- data.frame(Y = Y, f = factor_scores, Z = X[,2], X.1 = X[,3],
                           X.2 = X[,4], X.3 = X[,5])

  # Original outcome model plus the standardized orthogonal term
  outcome_model <- lm(Y ~ f * Z + I(f^2) * Z + X.1 + X.2 + X.3 + I(f^3) +
                        I(f^4), data = outcome_df)

  # Step 3: Generate data from estimated parameters
  N_full         <- nrow(X)
  N <- floor(N_factor * N_full)

  idx       <- sample(seq_len(N_full), size = N, replace = TRUE)
  new_X     <- X[idx,]
  new_X[,2] <- sample(0:1, size = N, replace = TRUE)
  bootstrapped_residuals <- X2_squared_orthogonal_std[idx]

  new_df <- data.frame(Z = new_X[,2], X.1 = new_X[,3],
                       X.2 = new_X[,4], X.3 = new_X[,5])

  # Generate latent factor based on covariates
  theta_mean <- predict(theta_model, new_df)
  theta_new <- theta_mean + sigma(theta_model) * rnorm(N)

  # Generate items
  M_sim <- matrix(NA, N, ncol(M))
  for(p in 1:ncol(M)) {
    M_sim[,p] <- item_params[[p]]$intercept +
      item_params[[p]]$loading * theta_new +
      rnorm(N, 0, item_params[[p]]$residual_sd)
  }

  new_df$f <- theta_new
  new_df$M <- M_sim

  # When predicting the outcome, add the bootstrapped nonlinear effect scaled by the strength parameter
  bootstrapped_residuals <- new_X[,4] * new_X[,5]
  Y_mean <- predict(outcome_model, new_df) +
    nonlinear_strength * bootstrapped_residuals

  Y_sim <- Y_mean + sigma_factor * sigma(outcome_model) * rnorm(N)
  new_df$Y <- Y_sim

  return(new_df)
}

# Function for composite score-based mediation analysis
# Uses factor analysis with regression scores
#' @export
fit_composite_mediation <- function(M, Y, Z, X) {
  # Step 1: Create composite score using factor analysis with regression scores
  fa_result <- factanal(M, factors = 1, scores = "regression")
  factor_scores <- fa_result$scores[,1]

    med_data <- data.frame(
      Y = Y,
      M = factor_scores,
      Z = Z,
      X
    )

    fit_M <- lm(M ~ . - Y, data = med_data)
    fit_Y <- lm(Y ~ ., data = med_data)

    mediation_results <- mediation::mediate(fit_M, fit_Y,
                                            mediator = "M", treat = "Z")


    return(mediation_results)
}
