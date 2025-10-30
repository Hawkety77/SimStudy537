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
generate_compound_symmetric_covariance_matrix(sigma2 = 10, rho = 0.5, n_measures = 5)

# Generate the random coefficient (X needs to be predetermined)
generate_random_coefficients_covariance_matrix(sigma_b2 = 1, sigma_c2 = 0.5, sigma_e2 = 10, x = matrix(c(0, 1, 2, 3, 4), nrow = 5, ncol = 1), n_measures = 5)

# Generate autoregressive covariance matrix
generate_autoregressive_covariance_matrix(sigma2 = 10, phi = 0.8, n_measures = 5)

