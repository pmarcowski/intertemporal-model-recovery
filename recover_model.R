# Title: Model Selection and Recovery
# Author: Przemyslaw Marcowski, PhD
# Email: p.marcowski@gmail.com
# Date: 2023-03-07
# Copyright (c) 2023 Przemyslaw Marcowski

# This code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

# Load packages
library(tidyverse)
library(ggridges)
library(see)
library(knitr)
library(patchwork)

# Local source
source("R/models.R")
source("R/tools.R")

# Set global parameters
set.seed(42) # for reproducibility
n_repetitions <- 100 # number of MCCV repetitions
n_fits <- 10 # number of refits per MCCV repetition
train_proportion <- 0.7 # proportion of train/test data split

# Checkpoint configuration
checkpoint_files <- list(
  cv_results = "checkpoints/cross_validation.rds",
  agents = "checkpoints/agent_parameters.rds",
  simulations = "checkpoints/simulations.rds",
  recovered = "checkpoints/recovered_parameters.rds"
)

# Data preparation --------------------------------------------------------
choices_data <- read_csv("data/choices.csv") %>%
  filter(Condition == 1) %>%
  group_by(Subject) %>%
  filter(all(!is.na(across(everything())))) %>%
  ungroup() %>%
  mutate(
    X_star = (X1 + X2) / 2,
    T_star = (T1 + T2) / 2,
    G = c(scale(X2 - X_star)),
    R = c(scale((X2 - X1) / X_star)),
    D = c(scale(T2 - T_star)),
    T = c(scale((T2 - T1) / T_star)),
    Drift_D = c(scale(X2 - X1)),
    Drift_R = c(scale((X2 - X1) / X1)),
    Drift_I = c(scale((X2 / X1)^(1 / (T2 - T1)) - 1)),
    Drift_T = c(scale(T2 - T1)),
    X1 = X1 / max(X2),
    X2 = X2 / max(X2)
  )

choices_list <- split(choices_data, choices_data$Subject, drop = TRUE)

# Model selection ---------------------------------------------------------

# Setup models
model_names <- c("EXPO", "HYPER2", "DEXPO", "ITCH", "DRIFT", "TRADE")
parameter_limits <- map(paste0(model_names, "_lim"), get)
names(parameter_limits) <- model_names
n_parameters <- sapply(model_names, function(m) ncol(parameter_limits[[m]]))

# Load or compute cross-validation results
cv_results <- load_checkpoint(checkpoint_files$cv_results)

if (is.null(cv_results)) {
  cv_results_list <- list()

  for (model in model_names) {
    start_time <- Sys.time()
    cat("\nEvaluate Model:", model, "- Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

    cat("Initializing parameters...")
    start_values <- get_starts(model, n_parameters[model], choices_list, parameter_limits[[model]])
    cat("DONE!", "\n")

    cat("Cross-validating...")
    results <- map(
      names(choices_list),
      ~ get_validate(
        ., choices_list, model, n_parameters[model],
        n_repetitions, n_fits, train_proportion, start_values,
        parameter_limits[[model]]
      )
    )

    cv_results_list[[model]] <- as.data.frame(do.call(rbind, results))
    cat("DONE!", "\n")

    end_time <- Sys.time()
    elapsed_time <- difftime(end_time, start_time, units = "secs")
    cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Elapsed:", round(elapsed_time, 3), "seconds", "\n")
  }

  # Combine validation results
  cv_results <- bind_rows(cv_results_list, .id = "model") %>%
    set_names(c(
      "model", "id", "n_fit", "n_test", "estim",
      "logloss", "mae", "extr", "n_na", "K",
      paste0("par", 1:max(n_parameters))
    )) %>%
    type.convert(as.is = TRUE)

  save_checkpoint(cv_results, checkpoint_files$cv_results)
}

# Performance summary
cv_results_summary <- cv_results %>%
  group_by(model) %>%
  summarise(
    mean_logloss = mean(logloss),
    pos = 0.1 + (quantile(logloss, 0.75) + 1.5 * IQR(logloss))
  ) %>%
  ungroup() %>%
  arrange(mean_logloss)

kable(
  cv_results_summary[-3],
  caption = "Cross-Validation Results",
  digits = 3
)

best_model <- cv_results_summary$model[1]

write_csv(cv_results_summary, "output/mean_model_error.csv")

# Create model error plot
model_error_plot <- cv_results %>%
  ggplot(aes(x = reorder(model, logloss), y = logloss, color = model, fill = model)) +
  geom_hline(yintercept = -log(0.5), linetype = "dashed") +
  geom_jitter(height = 0, width = 0.1, size = 0.5, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.2, color = "black") +
  stat_boxplot(geom = "errorbar", color = "black", width = 0.1) +
  stat_summary(fun = mean, geom = "point", shape = 1, size = 2, color = "black") +
  geom_text(
    data = cv_results_summary,
    aes(x = reorder(model, mean_logloss), y = pos, label = gsub("0\\.", "\\.", format(round(mean_logloss, 3)))),
    position = position_dodge(width = 0.75),
    color = "black"
  ) +
  labs(
    title = "Model Error",
    x = NULL, y = "Log-Loss",
    color = "Model", fill = "Model"
  ) +
  scale_color_see() +
  scale_fill_see() +
  theme_modern() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  annotate(
    "rect",
    xmin = which(levels(reorder(cv_results$model, cv_results$logloss)) == best_model) - 0.3,
    xmax = which(levels(reorder(cv_results$model, cv_results$logloss)) == best_model) + 0.3,
    ymin = -0.3, ymax = 3,
    color = "transparent", fill = "darkgreen", alpha = 0.2
  )

print(model_error_plot)
ggsave("output/model_error.png", model_error_plot)

# Model recovery ----------------------------------------------------------

# Load or compute agent parameters
agents <- load_checkpoint(checkpoint_files$agents)

if (is.null(agents)) {
  starts_agents <- get_starts(
    best_model, n_parameters[best_model],
    choices_list, parameter_limits[[best_model]]
  )

  agents <- map(choices_list, ~ {
    fit <- get_fit_optim(
      best_model, n_parameters[best_model],
      ., n_fits, starts_agents, parameter_limits[[best_model]]
    )

    fit[-1]
  })

  save_checkpoint(agents, checkpoint_files$agents)
}

# Load or compute simulations
simulations <- load_checkpoint(checkpoint_files$simulations)

if (is.null(simulations)) {
  # Simulate agents
  n_trials <- max(choices_data$Question)
  x_values <- sort(unique(c(choices_data$X1, choices_data$X2)))
  t_values <- sort(unique(c(choices_data$T1, choices_data$T2)))

  simulations_list <- map2(choices_list, agents, function(agent_data, agent_param) {
    simulate_agent_mel(agent_data, agent_param, get(best_model))
  })

  simulations <- bind_rows(simulations_list, .id = "id")

  save_checkpoint(simulations, checkpoint_files$simulations)
}

# Compare empirical and simulated choice distributions
empirical_choices <- choices_data$LaterOptionChosen
simulated_choices <- simulations$LaterOptionChosen

# Prepare data for plotting
choice_distribution <- bind_rows(
  tibble(
    group = "Empirical",
    choice = factor(empirical_choices, labels = c("Earlier", "Later"))
  ),
  tibble(
    group = "Simulated",
    choice = factor(simulated_choices, labels = c("Earlier", "Later"))
  )
) %>%
  group_by(group, choice) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup()

# Calculate rounded up maximum count
rounded_max_count <- ceiling(max(choice_distribution$count) / 1000) * 1000

# Create choice distribution plot
choice_plot <- choice_distribution %>%
  ggplot(aes(x = choice, y = count, fill = group)) +
  geom_bar(
    stat = "identity",
    width = 0.5, position = position_dodge(width = 0.5),
    color = "black", alpha = 0.8
  ) +
  geom_line(
    aes(y = proportion * rounded_max_count, group = group),
    position = position_dodge(width = 0.5), color = "black",
    linetype = "dashed"
  ) +
  geom_point(
    aes(y = proportion * rounded_max_count),
    position = position_dodge(width = 0.5)
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.5),
    vjust = -1
  ) +
  geom_text(
    aes(y = proportion * rounded_max_count, label = scales::percent(proportion, accuracy = 0.1)),
    position = position_dodge(width = 0.5),
    vjust = 2
  ) +
  scale_y_continuous(
    name = "Count",
    limits = c(0, rounded_max_count),
    breaks = seq(0, rounded_max_count, by = rounded_max_count / 5),
    sec.axis = sec_axis(~ . / rounded_max_count, name = "Proportion", labels = scales::percent)
  ) +
  scale_fill_see() +
  scale_color_see() +
  labs(
    title = "Choice Distribution Comparison",
    x = "Choice",
    fill = "Data Origin"
  ) +
  theme_modern() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

print(choice_plot)

# Create a 2x2 table of choices
choice_table <- table(Empirical = empirical_choices, Simulated = simulated_choices)
colnames(choice_table) <- c("Earlier Chosen", "Later Chosen")
rownames(choice_table) <- c("Empirical", "Simulated")

# Perform chi-square test
choice_test <- chisq.test(choice_table)
kable(choice_table)
print(choice_test)

ggsave("output/choice_distributions.png", choice_plot)
write_csv(as.data.frame(choice_table), "output/empirical_simulated_choices.csv")

# Load or compute recovered parameters
recovered <- load_checkpoint(checkpoint_files$recovered)

if (is.null(recovered)) {
  # Recover parameters by fitting model to simulated choices
  simulations_list <- split(simulations, simulations$id, drop = TRUE)

  starts_recovered <- get_starts(
    best_model, n_parameters[best_model],
    simulations_list, parameter_limits[[best_model]]
  )

  recovered <- map(simulations_list, ~ {
    fit <- get_fit_optim(
      best_model, n_parameters[best_model],
      ., n_fits, starts_recovered, parameter_limits[[best_model]]
    )

    fit[-1]
  })

  save_checkpoint(recovered, checkpoint_files$recovered)
}

# Prepare data for recovery analysis
empirical_parameters <- prepare_recovery_data(agents, "empirical", n_parameters[best_model])
recovered_parameters <- prepare_recovery_data(recovered, "recovered", n_parameters[best_model])

# Identify parameter columns (excluding 'id' and 'origin')
parameter_names <- setdiff(names(empirical_parameters), c("id", "origin"))

# Prepare data for parameter distribution comparison
parameter_distribution <- bind_rows(
  empirical_parameters %>% mutate(group = "Empirical"),
  recovered_parameters %>% mutate(group = "Recovered")
) %>%
  pivot_longer(cols = all_of(parameter_names), names_to = "parameter", values_to = "value")

# Create parameter distribution plot
parameter_plot <- parameter_distribution %>%
  ggplot(aes(x = value, y = parameter, fill = group)) +
  geom_density_ridges(
    alpha = 0.8, scale = 0.9,
    quantile_lines = TRUE, quantiles = 2,
    color = "black", linetype = "dashed"
  ) +
  scale_fill_see() +
  labs(
    title = "Parameter Distribution Comparison",
    x = "Value", y = "Parameter",
    fill = "Parameter Origin"
  ) +
  theme_modern() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

print(parameter_plot)
ggsave("output/parameter_distributions.png", parameter_plot)

# Calculate Bland-Altman values for each parameter
bland_altman_data <- map_dfr(parameter_names, ~ {
  x <- empirical_parameters[[.x]]
  y <- recovered_parameters[[.x]]

  mean_diff <- mean(x - y)
  sd_diff <- sd(x - y)

  tibble(
    par = .x,
    mean_value = (x + y) / 2,
    diff = x - y,
    mean_diff = mean_diff,
    upper_limit = mean_diff + 1.96 * sd_diff,
    lower_limit = mean_diff - 1.96 * sd_diff
  )
})

write_csv(bland_altman_data, "output/bland_altman_results.csv")

# Create Bland-Altman plot
bland_altman_plot <- bland_altman_data %>%
  ggplot(aes(x = mean_value, y = diff)) +
  geom_point(alpha = 0.5) +
  geom_hline(aes(yintercept = mean_diff), color = "black", linetype = "dashed") +
  geom_hline(aes(yintercept = upper_limit), color = "red3", linetype = "dashed") +
  geom_hline(aes(yintercept = lower_limit), color = "red3", linetype = "dashed") +
  labs(
    title = "Log-Transformed Bland-Altman Plots for Recovery Analysis",
    x = "Mean of Empirical and Recovered Parameters",
    y = "Parameter Difference\n(Empirical - Recovered)"
  ) +
  facet_wrap(~par, scales = "free", labeller = label_parsed) +
  theme_modern() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", fill = "transparent", linewidth = 0.5),
    strip.background = element_rect(color = "transparent", fill = "transparent")
  )

print(bland_altman_plot)
ggsave("output/parameter_recovery.png", bland_altman_plot, height = 5)

# Initialize results dataframe
results <- data.frame(
  parameter = parameter_names,
  observed_diff = numeric(length(parameter_names)),
  p_value = numeric(length(parameter_names))
)

# Perform permutation tests to compare empirical and recovered parameters
for (i in seq_along(parameter_names)) {
  param <- parameter_names[i]
  emp <- empirical_parameters[[param]]
  rec <- recovered_parameters[[param]]

  perm_test <- permutation_test(emp, rec)
  results$observed_diff[i] <- perm_test$observed_diff
  results$p_value[i] <- perm_test$p_value
}

# Inspect comparison results
kable(
  results,
  caption = "Comparison of Empirical and Recovered Parameters",
  digits = 3
)

write_csv(results, "output/empirical_simulated_parameters.csv")

# Figures -----------------------------------------------------------------

# Create figure for model evaluation and recovery
fig <- model_error_plot + parameter_plot + bland_altman_plot + choice_plot +
  plot_layout(widths = c(1, 0.8), heights = c(1, 0.8))

print(fig)
ggsave("output/figure.png", fig, height = 8, scale = 1.2)
