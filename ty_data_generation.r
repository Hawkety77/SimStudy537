# libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(nlme)

## NOTES
## id is the subject
## group is the treatment
## time is the repeated measure (i.e. what beta 1 is multiplied by)

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

# Example Generate Data
beta_0 <- 10
beta_1 <- 0
n_subjects_per_treatment <- 100

df_comp <- generate_data(beta_0, beta_1, n_subjects_per_treatment, V_comp)
df_random <- generate_data(beta_0, beta_1, n_subjects_per_treatment, V_random)
df_ar <- generate_data(beta_0, beta_1, n_subjects_per_treatment, V_ar)

# Example Compound Symetric
mod_full <- gls(y ~ 1 + time + time:group, data = df_comp,
                correlation = corCompSymm(form = ~ 1 | id),
                method = "ML")
p_lrt <- anova(mod_full)$"p-value"[3]

mod_reml <- gls(y ~ 1 + time + time:group, data = df_comp,
                correlation = corCompSymm(form = ~ 1 | id),
                method = "REML")
p_f <- anova(mod_reml)$"p-value"[3]
anova(mod_full)
anova(mod_reml)
print(p_lrt)
print(p_f)

# Example RC
# Example Autoregressive
mod_full <- gls(y ~ 1 + time + time:group, data = df_ar,
                correlation = corAR1(form = ~ time | id),
                method = "ML")
p_lrt <- anova(mod_full)$"p-value"[3]

mod_reml <- gls(y ~ 1 + time + time:group, data = df_ar,
                correlation = corAR1(form = ~ time | id),
                method = "REML")
p_f <- anova(mod_reml)$"p-value"[3]
anova(mod_full)
anova(mod_reml)
print(p_lrt)
print(p_f)

### Simulation ###
set.seed(77)
 
n_simulations_per_pval <- 100

beta_0 <- 10
beta_1 <- 1
alpha = .05
n_subjects_per_treatment_values <- c(3, 5, 10, 50)
covariance_structures <- list(V_comp, V_random, V_ar)
covariance_structures_names <- c("Compound Symmetric", "Random Coefficients", "Autoregressive")
results <- data.frame()
for (i in 1:n_simulations_per_pval){
  for (beta_1_sim in c(0, beta_1)){
    for (n in n_subjects_per_treatment_values){
      for (j in 1:3){
        print(paste("Simulation:", i, "Subjects:", n, "Covariance Structure:", covariance_structures_names[j], "Effect Size:", ifelse(beta_1_sim == 0, "null", "alternative")))
        df <- generate_data(beta_0, beta_1_sim, n, covariance_structures[[j]]) 
        # Fit models
        if (j == 1){
          mod_lrt <- gls(y ~ 1 + time + time:group, data = df,
                          correlation = corCompSymm(form = ~ 1 | id),
                          method = "ML")
          p_lrt <- anova(mod_lrt)$"p-value"[3]

          mod_reml <- gls(y ~ 1 + time + time:group, data = df,
                          correlation = corCompSymm(form = ~ 1 | id),
                          method = "REML")
          p_f <- anova(mod_reml)$"p-value"[3]
        }

        if (j == 2){
          mod_lrt <- lmer(y ~ 1 + time + time:group + (1 | id), data = df, REML = FALSE)
          p_lrt <- anova(mod_lrt)["time:group", "Pr(>F)"]

          mod_reml <- lmer(y ~ 1 + time + time:group + (1 | id), data=df, REML=TRUE)
          p_f <- anova(mod_reml)["time:group", "Pr(>F)"]
        }

        if (j == 3){
          mod_lrt <- gls(y ~ 1 + time + time:group, data = df,
                          correlation = corAR1(form = ~ time | id),
                          method = "ML")
          p_lrt <- anova(mod_lrt)$"p-value"[3]

          mod_reml <- gls(y ~ 1 + time + time:group, data = df,
                          correlation = corAR1(form = ~ time | id),
                          method = "REML")
          p_f <- anova(mod_reml)$"p-value"[3]
        }

        # Store results
        results <- rbind(results, data.frame(
          simulation = i,
          n_subjects = n,
          covariance_structure = covariance_structures_names[j],
          effect_size = ifelse(beta_1_sim == 0, "null", "alternative"), 
          p_lrt = round(p_lrt, 5),
          p_f = round(p_f, 5), 
          reject_lrt = ifelse(p_lrt < alpha, 1, 0),
          reject_f = ifelse(p_f < alpha, 1, 0), 
          accuracy_lrt = ifelse((p_lrt < alpha & beta_1_sim == 1) | (p_lrt > alpha & beta_1_sim == 0), 1, 0), 
          accuracy_f = ifelse((p_f < alpha & beta_1_sim == 1) | (p_f > alpha & beta_1_sim == 0), 1, 0)
        ))
      }
    }
  }
}

summary_results <- results %>%
  group_by(n_subjects, covariance_structure, effect_size) %>%
  summarize(
    power_LRT = mean(reject_lrt),
    power_F   = mean(reject_f),
    .groups = "drop"
  )

summary_long <- summary_results %>%
  pivot_longer(cols = c(power_LRT, power_F),
               names_to = "test",
               values_to = "prob") %>%
  mutate(effect_size = ifelse(effect_size == "null", "Type I Error", "Power"),
         test = ifelse(test == "power_LRT", "LRT", "F-test"))

summary_long_type1 <- summary_long %>% 
  filter(effect_size == "Type I Error")

ggplot(summary_long_type1,
       aes(x = n_subjects, y = prob,
           color = covariance_structure,
           linetype = test)) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal(base_size = 14) +
  labs(x = "Number of Subjects per Treatment",
       y = "Type I Error Rate",
       color = "Covariance Structure",
       linetype = "Test")
