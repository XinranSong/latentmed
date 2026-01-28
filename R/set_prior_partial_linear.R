set_prior_partial_linear <- function(
    J_M,       # Number of mediator items
    J_Y,       # Number of outcome items
    mu_gamma = 0,
    sigma_gamma2 = 5,
    mu_beta_M0 = NULL,
    sigma_beta_M02 = NULL,
    mu_beta_M1 = NULL,
    sigma_beta_M12 = NULL,
    mu_beta_Y0 = NULL,
    sigma_beta_Y02 = NULL,
    mu_beta_Y1 = NULL,
    sigma_beta_Y12 = NULL,
    a_Mp = 0.1,
    b_Mp = 0.1,
    a_Yq = 0.1,
    b_Yq = 0.1,
    a_M = 0.1,
    b_M = 0.1,
    a_Y = 0.1,
    b_Y = 0.1
) {
  # Set default priors for beta_M0 and beta_M1 if not provided
  ## if (is.null(mu_beta_M0))     mu_beta_M0     <- rep(0, J_M)
  ## if (is.null(sigma_beta_M02)) sigma_beta_M02 <- rep(100, J_M)
  ## if (is.null(mu_beta_M1))     mu_beta_M1     <- rep(1, J_M)
  ## if (is.null(sigma_beta_M12)) sigma_beta_M12 <- rep(100, J_M)
  if (is.null(mu_beta_M0))     mu_beta_M0     <- 0
  if (is.null(sigma_beta_M02)) sigma_beta_M02 <- 100
  if (is.null(mu_beta_M1))     mu_beta_M1     <- 1
  if (is.null(sigma_beta_M12)) sigma_beta_M12 <- 100

  # Set default priors for beta_Y0 and beta_Y1 if not provided
  ## if (is.null(mu_beta_Y0))     mu_beta_Y0     <- rep(0, J_Y)
  ## if (is.null(sigma_beta_Y02)) sigma_beta_Y02 <- rep(100, J_Y)
  ## if (is.null(mu_beta_Y1))     mu_beta_Y1     <- rep(1, J_Y)
  ## if (is.null(sigma_beta_Y12)) sigma_beta_Y12 <- rep(100, J_Y)
  if (is.null(mu_beta_Y0))     mu_beta_Y0     <- 0
  if (is.null(sigma_beta_Y02)) sigma_beta_Y02 <- 100
  if (is.null(mu_beta_Y1))     mu_beta_Y1     <- 1
  if (is.null(sigma_beta_Y12)) sigma_beta_Y12 <- 100

  # Check that the lengths of the vectors match the number of items
  ## if(length(mu_beta_M0) != J_M) stop("Length of mu_beta_M0 must be equal to J_M")
  ## if(length(sigma_beta_M02) != J_M) stop("Length of sigma_beta_M02 must be equal to J_M")
  ## if(length(mu_beta_M1) != J_M) stop("Length of mu_beta_M1 must be equal to J_M")
  ## if(length(sigma_beta_M12) != J_M) stop("Length of sigma_beta_M12 must be equal to J_M")
  ## if(length(mu_beta_Y0) != J_Y) stop("Length of mu_beta_Y0 must be equal to J_Y")
  ## if(length(sigma_beta_Y02) != J_Y) stop("Length of sigma_beta_Y02 must be equal to J_Y")
  ## if(length(mu_beta_Y1) != J_Y) stop("Length of mu_beta_Y1 must be equal to J_Y")
  ## if(length(sigma_beta_Y12) != J_Y) stop("Length of sigma_beta_Y12 must be equal to J_Y")

  priors <- list(
    # Prior for gamma
    mu_gamma = mu_gamma,             # Prior mean for gamma
    sigma_gamma2 = sigma_gamma2,     # Prior variance for gamma

    # Priors for beta_M parameters
    mu_beta_M0 = mu_beta_M0,         # Prior means for beta^M_{0j}
    sigma_beta_M02 = sigma_beta_M02, # Prior variances for beta^M_{0j}
    mu_beta_M1 = mu_beta_M1,         # Prior means for beta^M_{1j}
    sigma_beta_M12 = sigma_beta_M12, # Prior variances for beta^M_{1j}

    # Priors for beta_Y parameters
    mu_beta_Y0 = mu_beta_Y0,         # Prior means for beta^Y_{0j}
    sigma_beta_Y02 = sigma_beta_Y02, # Prior variances for beta^Y_{0j}
    mu_beta_Y1 = mu_beta_Y1,         # Prior means for beta^Y_{1j}
    sigma_beta_Y12 = sigma_beta_Y12, # Prior variances for beta^Y_{1j}

    # Priors for variance parameters (inverse-gamma hyperparameters)
    a_Mp = a_Mp,  # Shape parameter for a^2_j (mediator measurement error variances)
    b_Mp = b_Mp,  # Scale parameter for a^2_j
    a_Yq = a_Yq,  # Shape parameter for b^2_j (outcome measurement error variances)
    b_Yq = b_Yq,  # Scale parameter for b^2_j
    a_M = a_M,    # Shape parameter for sigma^2_M (variance of theta_M)
    b_M = b_M,    # Scale parameter for sigma^2_M
    a_Y = a_Y,    # Shape parameter for sigma^2_Y (variance of theta_Y)
    b_Y = b_Y     # Scale parameter for sigma^2_Y
  )

  return(priors)
}
