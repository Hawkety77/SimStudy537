library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)

results <- read.csv("simulation_results.csv")

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
