# libraries
library(tidyverse)

generate_compound_symmetric_covariance_matrix <- function(sigma2, rho, n_measures){
  V_within <- sigma2 * ((1 - rho) * diag(n_measures) + rho * matrix(1, n_measures, n_measures))
  return(V_within)
}

generate_random_coefficients_covariance_matrix <- function(sigma_b2, sigma_c2, sigma_e2, x, n_measures){
  ones <- matrix(1, nrow = n_measures, ncol = 1)
  I <- diag(n_measures)
  V <- sigma_b2 * (ones %*% t(ones)) + sigma_c2 * (x %*% t(x)) + sigma_e2 * I
  return(V)
}

generate_autoregressive_covariance_matrix <- function(sigma2, phi, n_measures){
  V <- matrix(0, n_measures, n_measures)
  for(i in 1:n_measures){
    for(j in 1:n_measures){
      V[i, j] <- sigma2 * phi^abs(i - j)
    }
  }
  return(V)
}

# Generate the compound symmetric covariance matrix
V_comp <- generate_compound_symmetric_covariance_matrix(sigma2 = 10, rho = 0.5, n_measures = 5)

# Generate the random coefficient ***(X is predetermined here!!!)***
V_random <- generate_random_coefficients_covariance_matrix(sigma_b2 = 1, sigma_c2 = 0.5, sigma_e2 = 5, x = matrix(c(0, 1, 2, 3, 4), nrow = 5, ncol = 1), n_measures = 5)

# Generate autoregressive covariance matrix
V_ar <- generate_autoregressive_covariance_matrix(sigma2 = 10, phi = 0.584776, n_measures = 5)

# Check that the determinants are approximately the same
print(det(V_comp))
print(det(V_random))
print(det(V_ar))

# Generate data function
generate_data <- function(beta_0, beta_1, n_subjects_per_treatment, V){

  # Same between treatments
  V_full <- kronecker(diag(n_subjects_per_treatment), V)
  C <- chol(V_full)
  X <- cbind(1, rep(1:5, times = n_subjects_per_treatment))

  # Control
  Z <- matrix(rep(rnorm(n_subjects_per_treatment, 0, 1), each = 5), ncol = 1, byrow = TRUE)
  beta <- matrix(c(beta_0, 0), nrow = 2, byrow = TRUE)
  y_control <- X %*% beta + C %*% Z

  # Treatment
  Z <- matrix(rep(rnorm(n_subjects_per_treatment, 0, 1), each = 5), ncol = 1, byrow = TRUE)
  beta <- matrix(c(beta_0, beta_1), nrow = 2, byrow = TRUE)
  y_treatment <- X %*% beta + C %*% Z

  data <- data.frame(
    id = rep(1:(n_subjects_per_treatment* 2), each = 5),
    time = rep(1:5, times = n_subjects_per_treatment),
    group = rep(c("0", "1"), each = 5 * n_subjects_per_treatment), 
    y = c(y_control, y_treatment)
  )

  return(data)

}

# Generate Data
beta_0 <- 10
beta_1 <- 5
n_subjects_per_treatment <- 4

generate_data(beta_0, beta_1, n_subjects_per_treatment, V_comp)
generate_data(beta_0, beta_1, n_subjects_per_treatment, V_random)
generate_data(beta_0, beta_1, n_subjects_per_treatment, V_ar)



