# libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(nlme)

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

V_comp <- generate_compound_symmetric_covariance_matrix(sigma2 = 10, rho = 0.5, n_measures = 5)
V_random <- generate_random_coefficients_covariance_matrix(sigma_b2 = 1, sigma_c2 = 1, sigma_e2 = 4.3833424, x = matrix(c(0, 1, 2, 3, 4), nrow = 5, ncol = 1), n_measures = 5)
V_ar <- generate_autoregressive_covariance_matrix(sigma2 = 10, phi = 0.584776, n_measures = 5)

# Check that the determinants are approximately the same
print(det(V_comp))
print(det(V_random))
print(det(V_ar))

# Generate data function
generate_data <- function(beta_0, beta_1, n_subjects_per_treatment, V) {
  n_subjects <- n_subjects_per_treatment * 2
  n_time <- ncol(V)
  variance <- MASS::mvrnorm(n = n_subjects, mu = rep(0, n_time), Sigma = V)
  df <- as.data.frame(variance) %>%
    mutate(id = 1:n_subjects) %>%
    pivot_longer(-id, names_to = "time", values_to = "error") |>
    mutate(
      time = rep(0:4, times = n_subjects) ,
      treatment = ifelse(id <= n_subjects_per_treatment, 0, 1),
      y = beta_0 + beta_1 * time * treatment + error
    ) %>%
    select(id, treatment, time, y)
  
  return(df)
}

df_random <- generate_data(10, 5, 10, V_random)

ggplot(df_random, aes(x = time, y = y, color = treatment, group = id)) +
  geom_line(alpha = 0.3) +
  theme_minimal()


### Simulation ###
set.seed(77)
 
n_simulations_per_pval <- 10000

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
          mod_lrt <- gls(y ~ 1 + time + time:treatment, data = df,
                          correlation = corCompSymm(form = ~ 1 | id),
                          method = "ML")
          p_lrt <- anova(mod_lrt)$"p-value"[3]

          mod_reml <- gls(y ~ 1 + time + time:treatment, data = df,
                          correlation = corCompSymm(form = ~ 1 | id),
                          method = "REML")
          p_f <- anova(mod_reml)$"p-value"[3]
        }

        if (j == 2){
          mod_lrt <- lmer(y ~ 1 + time + time:treatment + (1 + time | id), data=df, REML=FALSE)
          p_lrt <- anova(mod_lrt, type = 3)[["Pr(>F)"]][2]

          mod_reml <- lmer(y ~ 1 + time + time:treatment + (1 + time | id), data=df, REML=TRUE)
          p_f <- anova(mod_reml, type = 3)[["Pr(>F)"]][2]
        }

        if (j == 3){
          mod_lrt <- gls(y ~ 1 + time + time:treatment, data = df,
                          correlation = corAR1(form = ~ time | id),
                          method = "ML")
          p_lrt <- anova(mod_lrt)$"p-value"[3]

          mod_reml <- gls(y ~ 1 + time + time:treatment, data = df,
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

# save results to a .csv
write.csv(results, "simulation_results.csv", row.names = FALSE)

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

# do power now
summary_long_power <- summary_long %>% 
  filter(effect_size == "Power")
ggplot(summary_long_power,
       aes(x = n_subjects, y = prob,
           color = covariance_structure,
           linetype = test)) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal(base_size = 14) +
  labs(x = "Number of Subjects per Treatment",
       y = "Power",
       color = "Covariance Structure",
       linetype = "Test")


# make table of results
library(knitr)
library(kableExtra)

# Type I Error table
type1_table <- summary_results %>%
  filter(effect_size == "null") %>%
  pivot_longer(cols = c(power_LRT, power_F),
               names_to = "test",
               values_to = "prob") %>%
  mutate(test = ifelse(test == "power_LRT", "LRT", "F-test")) %>%
  pivot_wider(names_from = c(covariance_structure, test), values_from = prob)

kable(type1_table, digits = 3, caption = "Type I Error Rates by Covariance Structure and Test") %>%
  kable_styling(full_width = FALSE)

# Power table
power_table <- summary_results %>%
  filter(effect_size == "alternative") %>%
  pivot_longer(cols = c(power_LRT, power_F),
               names_to = "test",
               values_to = "prob") %>%
  mutate(test = ifelse(test == "power_LRT", "LRT", "F-test")) %>%
  pivot_wider(names_from = c(covariance_structure, test), values_from = prob)

kable(power_table, digits = 3, caption = "Power by Covariance Structure and Test") %>%
  kable_styling(full_width = FALSE)
