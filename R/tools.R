# Title: Helper functions for model evaluation
# Author: Przemyslaw Marcowski, PhD
# Email: p.marcowski@gmail.com
# Date: 2023-05-12
# Copyright (c) 2023 Przemyslaw Marcowski

# This code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

# This script contains helper functions for model evaluation and simulation

save_checkpoint <- function(data, file_path) {
  # Saves data as an RDS file at the specified file path.
  #
  # Args:
  #   data: The data to be saved.
  #   file_path: The file path where the data should be saved.
  #
  # Returns:
  #   None. Saves the data and prints a message as a side effect.
  
  saveRDS(data, file_path)
  cat("Checkpoint saved:", file_path, "\n")
}

load_checkpoint <- function(file_path) {
  # Loads data from an RDS file if it exists, otherwise returns NULL.
  #
  # Args:
  #   file_path: The file path from which to load the data.
  #
  # Returns:
  #   The loaded data if the file exists, NULL otherwise.
  
  if (file.exists(file_path)) {
    cat("Loading checkpoint:", file_path, "\n")
    readRDS(file_path)
  } else {
    NULL
  }
}

# Helper function for preparing recovery analysis data
prepare_recovery_data <- function(parameters, origin, n_parameters) {
  # Prepares parameter data for recovery analysis by combining, renaming, and transforming.
  #
  # Args:
  #   parameters: A list of parameter sets.
  #   origin: A string indicating the origin of the data ("empirical" or "simulated").
  #   n_parameters: The number of parameters for the best model.
  #
  # Returns:
  #   A data frame with prepared parameter data, including an id column,
  #   renamed parameter columns, log-transformed parameter values, and a origin column.
  
  # Local function for offset and log-transformation
  .offset_log <- function(x) {
    # Applies offset and log-transformation to avoid negative values and stabilize variances.
    #
    # Args:
    #   x: Numeric vector to be transformed.
    #
    # Returns:
    #   A numeric vector with offset applied and log-transformed.
    #
    # Details:
    #   Adds a small offset (1e-7) to the absolute minimum value before log transformation
    #   to handle zero or negative values. This ensures all values are positive for logging.
    
    offset <- abs(min(x)) + 1e-7
    y <- x + offset
    return(log(y))
  }
  
  as.data.frame(do.call(rbind, parameters)) %>%
    mutate(id = names(parameters)) %>%
    select(id, everything()) %>%
    set_names(c("id", paste0("par", 1:n_parameters))) %>%
    mutate(across(-id, .offset_log)) %>%
    mutate(origin = origin)
}

# Data splitting function
split_data_mccv <- function(data, prop) {
  # Subsamples data into training and test set according to proportion prop.
  # Performs Monte Carlo Cross-Validation if prop < 1.
  #
  # Args:
  #   data: Data to split.
  #   prop: Proportion for splitting data into training and test set.
  #
  # Returns:
  #   A list containing training and test data sets. Test set is equal to training if prop = 1.
  
  if (prop < 1) {
    n_choices <- nrow(data)
    n_train <- floor(n_choices * prop)
    train_indices <- sample(n_choices, n_train)
    list(train = data[train_indices, ], test = data[-train_indices, ])
  } else {
    list(train = data, test = data)
  }
}

# Parameter initialization functions
draw_starts <- function(params) {
  # Draws start values for each model parameter from existing set.
  #
  # Args:
  #   params: Matrix of existing parameter sets.
  #
  # Returns:
  #   A numeric vector of starting values for each model parameter drawn from set.
  
  vapply(seq_len(ncol(params)), function(i) sample(params[, i], 1), numeric(1))
}

rand_starts <- function(model_name, param_limits) {
  # Randomizes start values for each model parameter within pre-defined limits.
  #
  # Args:
  #   model_name: Model string for which parameters are randomized.
  #   param_limits: List of parameter constraints for each model.
  #
  # Returns:
  #   A numeric vector of randomized starting values for each model parameter.
  
  vapply(seq_len(ncol(param_limits)), function(i) {
    runif(1,
          max(param_limits[1, i], -1e+10),
          min(param_limits[2, i], 1e+10)
    )
  }, numeric(1))
}

get_starts <- function(model_name, num_parameters, data_list, param_limits) {
  # Fits a model to data once to obtain a set of starting values.
  #
  # Args:
  #   model_name: Model string for which parameters are to be estimated.
  #   num_parameters: Number of model parameters.
  #   data_list: List of dataframes.
  #   param_limits: List of parameter constraints for each model.
  #
  # Returns:
  #   A matrix of parameter values.
  
  parameters <- matrix(NA, nrow = length(data_list), ncol = num_parameters)
  
  for (i in seq_along(data_list)) {
    fit <- FALSE
    while (!fit) {
      res_tmp <- tryCatch({
        optim(
          rand_starts(model_name, param_limits), get(model_name), dat = data_list[[i]],
          method = ifelse(num_parameters > 1, "L-BFGS-B", "Brent"),
          lower = param_limits[1, ], upper = param_limits[2, ]
        )
      }, error = function(e) {
        message("Error in optimization: ", e$message)
        return(NULL)
      })
      
      if (!is.null(res_tmp)) fit <- TRUE
    }
    parameters[i, ] <- res_tmp$par
  }
  
  parameters
}

# Model fitting function
get_fit_optim <- function(model_name, num_parameters, fit_data, num_fits, start_vals, param_limits) {
  # Fits model using L-BFGS-B (or Brent for one-dimensional models).
  # Fitting is performed by estimator MINIMIZATION.
  # Fitting is repeated until num_fits successful fits are obtained.
  #
  # Args:
  #   model_name: Model string for model to be fitted.
  #   num_parameters: Number of model parameters.
  #   fit_data: Data for fitting.
  #   num_fits: Number of fits performed to get best fit.
  #   start_vals: List of starting parameter values for each model.
  #   param_limits: List of parameter constraints for each model.
  #
  # Returns:
  #   A numeric vector containing the estimator value resulting from best fit,
  #   followed by the parameter values resulting from best fit.
  
  best_fit_value <- Inf
  successful_fits <- 0
  total_fits <- 0
  best_fit_params <- NULL
  
  while (successful_fits < num_fits && total_fits < 100) {
    total_fits <- total_fits + 1
    
    res_tmp <- tryCatch({
      optim(
        draw_starts(start_vals), get(model_name), dat = fit_data,
        method = ifelse(num_parameters > 1, "L-BFGS-B", "Brent"),
        lower = param_limits[1, ], upper = param_limits[2, ]
      )
    }, error = function(e) {
      message("Error in optimization: ", e$message)
      return(NULL)
    })
    
    if (!is.null(res_tmp)) {
      successful_fits <- successful_fits + 1
      if (res_tmp$value < best_fit_value) {
        best_fit_value <- res_tmp$value
        best_fit_params <- res_tmp$par
      }
    }
  }
  
  if (is.finite(best_fit_value)) {
    c(best_fit_value, best_fit_params)
  } else {
    NULL
  }
}

# Results extraction function
get_results_discrete <- function(test_data, model_name, fit) {
  # Extracts model training and test results in classification problems.
  #
  # Args:
  #   test_data: Test data set.
  #   model_name: Model string for model to be evaluated.
  #   fit: Model fit (vector containing estimator value and parameter values).
  #   epsilon: Smoothing parameter to avoid absolute certainty.
  #
  # Returns:
  #   A numeric vector containing:
  #   - Training model log-likelihood.
  #   - Multiple loss functions values obtained based on test data set.
  #   - Extremity abs(pred - .5).
  #   - Number of NAs obtained during modeling.
  
  estimator_value <- fit[1]
  pred <- do.call(model_name, list(fit[-1], test_data, choice = TRUE))
  
  logloss <- -mean(log(ifelse(test_data$LaterOptionChosen == 1, pred, 1 - pred)))
  mean_abs_error <- mean(abs(pred - test_data$LaterOptionChosen), na.rm = TRUE)
  extr <- mean(abs(pred - 0.5), na.rm = TRUE)
  n_na <- sum(is.na(pred))
  
  c(estimator_value, logloss, mean_abs_error, extr, n_na)
}

# Model validation function
get_validate <- function(id_name, data_list, model_name, num_parameters, num_repetitions, num_fits, split_prop, start_vals, param_limits) {
  # Evaluates a model for data set num_repetitions times.
  # Trains model using training data. Validates model using test data.
  # Models are evaluated based on the best of num_fits fit attempts.
  #
  # Args:
  #   id_name: Name or ID of the subject.
  #   data_list: List of dataframes containing all subjects' data.
  #   model_name: Model string for model to be evaluated.
  #   num_parameters: Number of model parameters.
  #   num_repetitions: Number of times to repeat the evaluation.
  #   num_fits: Number of fits performed to get best fit.
  #   split_prop: Proportion for splitting data into training and test set.
  #   start_vals: List of starting parameter values for each model.
  #   param_limits: List of parameter constraints for each model.
  #
  # Returns:
  #   A matrix containing:
  #   - Number of training and test set data points.
  #   - Training log-likelihoods.
  #   - Test set loss functions values for each repetition.
  #   - Number of estimated model parameters.
  #   - Parameter estimate values.
  
  metadata <- matrix(NA, ncol = 3, nrow = num_repetitions)
  results <- matrix(NA, ncol = 5, nrow = num_repetitions)
  parameters <- matrix(NA, ncol = num_parameters, nrow = num_repetitions)
  
  for (i in seq_len(num_repetitions)) {
    d_sets <- split_data_mccv(data_list[[id_name]], split_prop)
    metadata[i, ] <- c(id_name, nrow(d_sets$train), nrow(d_sets$test))
    fit <- get_fit_optim(model_name, num_parameters, d_sets$train, num_fits, start_vals, param_limits)
    results[i, ] <- get_results_discrete(d_sets$test, model_name, fit)
    parameters[i, ] <- fit[-1]
  }
  
  cbind(metadata, results, num_parameters, parameters)
}

# Agent simulation function
simulate_agent_mel <- function(agent_data, agent_param, choice_model) {
  # Simulates agent choices for MEL procedure.
  #
  # Args:
  #   n_trials: Number of choice trials to simulate.
  #   x_range: Numeric vector of min and max values for X1 and X2.
  #   t_range: Numeric vector of min and max values for T1 and T2.
  #   agent_param: Numeric vector of parameters for the agent.
  #   choice_model: Model function returning choice probabilities.
  #
  # Returns:
  #   A data frame containing simulated choices for the agent.
  
  results <- vapply(seq_len(nrow(agent_data)), function(trial) {
    X1 <- agent_data$X1[trial]
    X2 <- agent_data$X2[trial]
    T1 <- agent_data$T1[trial]
    T2 <- agent_data$T2[trial]
    
    pred <- do.call(
      choice_model,
      list(
        agent_param,
        list(X1 = X1, X2 = X2, T1 = T1, T2 = T2),
        choice = TRUE
      )
    )
    
    prob_later <- as.numeric(pred[1])
    choice_later <- rbinom(1, 1, prob_later)
    
    c(
      trial = trial,
      X1 = X1,
      X2 = X2,
      T1 = T1,
      T2 = T2,
      LaterOptionProb = prob_later,
      LaterOptionChosen = choice_later
    )
  }, numeric(7))
  
  as.data.frame(t(results))
}

# Permutation test function
permutation_test <- function(x, y, n_permutations = 10000) {
  # Performs a two-sided permutation test for paired samples.
  # It calculates the observed mean difference between x and y, then generates
  # a null distribution by randomly permuting the combined data n_permutations times.
  # The p-value is the proportion of permuted differences that are as extreme as
  # or more extreme than the observed difference in absolute value.
  #
  # Args:
  #   x: Numeric vector of the first sample.
  #   y: Numeric vector of the second sample.
  #   n_permutations: Number of permutations to perform (default: 10000).
  #
  # Returns:
  #   A list containing:
  #     observed_diff: The observed mean difference between x and y.
  #     p_value: The p-value from the permutation test.
  #
  
  observed_diff <- mean(x - y)
  
  combined <- c(x, y)
  n <- length(x)
  
  permuted_diffs <- replicate(n_permutations, {
    permuted <- sample(combined)
    mean(permuted[1:n] - permuted[(n + 1):(2 * n)])
  })
  
  p_value <- mean(abs(permuted_diffs) >= abs(observed_diff))
  
  list(observed_diff = observed_diff, p_value = p_value)
}
