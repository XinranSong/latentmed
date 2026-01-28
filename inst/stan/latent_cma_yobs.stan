data {
  int<lower=0> J;          // number of observations
  int<lower=0> P;          // number of items for mediator
  int<lower=0> Q;          // number of items for outcome
  int<lower=0> K;          // number of covariates
  matrix[J, K] X;          // covariates for latent variables
  matrix[J, P] M;          // observed mediator items
  vector[J] Y;             // observed outcome items
}

parameters {
  vector[K] alpha_M;        // coefficients for latent mediator model
  vector[K] alpha_Y;        // coefficients for outcome model
  real gamma;              // effect of latent mediator on outcome

  vector[J] theta_M;        // latent mediator variable

  vector[P] b_M0;           // uncentered intercepts for mediator items
  vector[P] b_M1;           // uncentered slopes for mediator items

  real<lower=0> sigma_M2;   // variance for latent mediator
  real<lower=0> sigma_Y2;   // variance for outcome
  vector<lower=0>[P] sigma_Mp2;  // item-specific variances for mediator items
}

transformed parameters {
  real<lower=0> sigma_Y;
  real<lower=0> sigma_M;
  real zeta;
  real delta;
  real tau;
  // Centering and normalizing intercepts and slopes
  vector[P] beta_M0;  // Centered intercepts for mediator items
  vector[P] beta_M1;  // Normalized slopes for mediator items

  // linking
  beta_M0 = b_M0;
  beta_M1 = b_M1;
  sigma_Y = sigma_Y2;
  sigma_M = sigma_M2;

  // Constraints for beta_M0 and beta_Y0
  // Applying identifiability constraints
  beta_M0[1] = 0;  // Constraint beta_M0[1] to be 0
  beta_M1[1] = 1;  // Constraint beta_M1[1] to be 1

   
  // Latent variable models
  vector[J] mu_M = X * alpha_M;  // mean of latent mediator
  vector[J] mu_Y = X * alpha_Y + gamma * theta_M;  // mean of latent outcome
  zeta = alpha_Y[2];
  delta = gamma * alpha_M[2];
  tau = zeta + delta;
}

model {
  // Priors for regression coefficients
  alpha_M ~ normal(0, 5);  // sigma_alpha_M2 = 1
  alpha_Y ~ normal(0, 5);  // sigma_alpha_Y2 = 1
  gamma ~ normal(0, 5);    // sigma_gamma2 = 1

  // Priors for item response parameters
  b_M0 ~ normal(0, 5);   // sigma_beta_M02 = 0.5
  b_M1 ~ normal(0, 5);   // sigma_beta_M12 = 0.5

  // Priors for variances
  sigma_M2 ~ uniform(0,10);
  sigma_Y2 ~ uniform(0, 10);
  sigma_Mp2 ~ uniform(0,10);

  // Likelihood for latent variables
  theta_M ~ normal(mu_M, sigma_M);
  Y ~ normal(mu_Y, sigma_Y);

  // Likelihood for observed mediator items
  for (j in 1:J) {
    for (p in 1:P) {
      M[j, p] ~ normal(beta_M0[p] + beta_M1[p] * theta_M[j], sigma_Mp2[p]);
    }
  }
}
