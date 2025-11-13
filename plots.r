library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)

results <- read.csv("SimStudy537/simulation_results.csv")

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
  geom_line(linewidth = .5) +
  scale_y_continuous(limits = c(0,.20)) +
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
  geom_line(linewidth = .5) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal(base_size = 14) +
  labs(x = "Number of Subjects per Treatment",
       y = "Power",
       color = "Covariance Structure",
       linetype = "Test")


# make table of results

# Type I Error Frequency Table
type1_table <- summary_results %>%
  filter(effect_size == "null") %>%
  select(n_subjects, covariance_structure, power_LRT, power_F) %>%
  pivot_longer(cols = c(power_LRT, power_F),
               names_to = "test",
               values_to = "prob") %>%
  mutate(test = ifelse(test == "power_LRT", "LRT", "REML")) %>%
  pivot_wider(names_from = c(covariance_structure, test),
              values_from = prob) %>%
  arrange(n_subjects)

colnames(type1_table) <- c("n_subjects", "LRT", "REML", "LRT", "REML", "LRT", "REML")

kable(type1_table, digits = 3, caption = "Type I Error Frequency") %>%
  add_header_above(c(" " = 1,
                     "Autoregressive" = 2,
                     "Random" = 2,
                     "Compound Symmetric" = 2)) %>%
  kable_styling(full_width = FALSE)

# Successful Rejection Frequency (Power) Table
power_table <- summary_results %>%
  filter(effect_size == "alternative") %>%
  select(n_subjects, covariance_structure, power_LRT, power_F) %>%
  pivot_longer(cols = c(power_LRT, power_F),
               names_to = "test",
               values_to = "prob") %>%
  mutate(test = ifelse(test == "power_LRT", "LRT", "REML")) %>%
  pivot_wider(names_from = c(covariance_structure, test),
              values_from = prob) %>%
  arrange(n_subjects)

colnames(power_table) <- c("n_subjects", "LRT", "REML", "LRT", "REML", "LRT", "REML")

kable(power_table, digits = 3, caption = "Successful Rejection Frequency") %>%
  add_header_above(c(" " = 1,
                     "Autoregressive" = 2,
                     "Random" = 2,
                     "Compound Symmetric" = 2)) %>%
  kable_styling(full_width = FALSE)

##

