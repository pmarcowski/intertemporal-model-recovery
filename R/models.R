# Title: Model functions for model evaluation
# Author: Przemyslaw Marcowski, PhD
# Email: p.marcowski@gmail.com
# Date: 2023-05-12
# Copyright (c) 2023 Przemyslaw Marcowski

# This code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

# This script contains model functions for model evaluation and simulation

EXPO_lim <- matrix(c(1e-3, 1e+2, -1e+1, 1e+1), nrow = 2)
EXPO <- function(par, dat, choice = FALSE, epsilon = 0.001, lambda = 0.001) {
  # Exponential model function
  #
  # Args:
  #   par: Numeric vector of parameters.
  #   dat: Data frame with necessary columns X1, X2, T1, and T2.
  #   choice: Boolean indicating whether to return choice probability.
  #   epsilon: Numeric smoothing parameter between 0 and 1.
  #   lambda: Numeric regularization parameter.
  #
  # Returns:
  #   If choice=TRUE, a numeric vector of choice probabilities.
  #   If choice=FALSE, a numeric representing the negative log-likelihood.
  
  # Compute option values
  ull <- dat$X2 * par[1]^dat$T2
  uss <- dat$X1 * par[1]^dat$T1
  
  # Compute choice probability
  P <- ull^par[2] / (ull^par[2] + uss^par[2])
  
  # Apply smoothing
  P <- epsilon * 0.5 + (1 - epsilon) * P
  
  # If choice=TRUE, return probabilities
  if (choice) {
    return(P)
  } else {
    # Otherwise, compute and return negative log posterior
    P <- ifelse(as.logical(dat$LaterOptionChosen), P, 1 - P)
    nll <- -sum(log(P))
    nlpr <- -sum(dnorm(par, 0, 10, log = TRUE)) # negative log prior
    reg <- lambda * sum(par^2) # regularization term
    nlpo <- nll + nlpr + reg # negative log posterior
    return(nlpo)
  }
}

HYPER2_lim <- matrix(c(1e-3, 1e+1, 1e-7, 1e+1, -1e+1, 1e+1), nrow = 2)
HYPER2 <- function(par, dat, choice = FALSE, epsilon = 0.001, lambda = 0.001) {
  # Hyperbolic model function
  #
  # Args:
  #   par: Numeric vector of parameters.
  #   dat: Data frame with necessary columns X1, X2, T1, and T2.
  #   choice: Boolean indicating whether to return choice probability.
  #   epsilon: Numeric smoothing parameter between 0 and 1.
  #
  # Returns:
  #   If choice=TRUE, a numeric vector of choice probabilities.
  #   If choice=FALSE, a numeric representing the negative log-likelihood.
  
  # Compute option values
  ull <- dat$X2 * (1 / (1 + par[1] * dat$T2)^par[2])
  uss <- dat$X1 * (1 / (1 + par[1] * dat$T1)^par[2])
  
  # Compute choice probability
  P <- ull^par[3] / (ull^par[3] + uss^par[3])
  
  # Apply smoothing
  P <- epsilon * 0.5 + (1 - epsilon) * P
  
  # If choice=TRUE, return probabilities
  if (choice) {
    return(P)
  } else {
    # Otherwise, compute and return negative log posterior
    P <- ifelse(as.logical(dat$LaterOptionChosen), P, 1 - P)
    nll <- -sum(log(P))
    nlpr <- -sum(dnorm(par, 0, 10, log = TRUE)) # negative log prior
    reg <- lambda * sum(par^2) # regularization term
    nlpo <- nll + nlpr + reg # negative log posterior
    return(nlpo)
  }
}

DEXPO_lim <- matrix(c(1e-3, 1, 1e-3, 1, 1e-7, 1, -1e+1, 1e+1), nrow = 2)
DEXPO <- function(par, dat, choice = FALSE, epsilon = 0.001, lambda = 0.001) {
  # Dual-Exponential model function
  #
  # Args:
  #   par: Numeric vector of parameters.
  #   dat: Data frame with necessary columns X1, X2, T1, and T2.
  #   choice: Boolean indicating whether to return choice probability.
  #   epsilon: Numeric smoothing parameter between 0 and 1.
  #
  # Returns:
  #   If choice=TRUE, a numeric vector of choice probabilities.
  #   If choice=FALSE, a numeric representing the negative log-likelihood.
  
  # Compute option values
  ull <- dat$X2 * ((par[3] * par[1]^dat$T2) + ((1 - par[3]) * par[2]^dat$T2))
  uss <- dat$X1 * ((par[3] * par[1]^dat$T1) + ((1 - par[3]) * par[2]^dat$T1))
  
  # Compute choice probability
  P <- ull^par[4] / (ull^par[4] + uss^par[4])
  
  # Apply smoothing
  P <- epsilon * 0.5 + (1 - epsilon) * P
  
  # If choice=TRUE, return probabilities
  if (choice) {
    return(P)
  } else {
    # Otherwise, compute and return negative log posterior
    P <- ifelse(as.logical(dat$LaterOptionChosen), P, 1 - P)
    nll <- -sum(log(P))
    nlpr <- -sum(dnorm(par, 0, 10, log = TRUE)) # negative log prior
    reg <- lambda * sum(par^2) # regularization term
    nlpo <- nll + nlpr + reg # negative log posterior
    return(nlpo)
  }
}

ITCH_lim <- matrix(c(-Inf, Inf, -Inf, Inf, -Inf, Inf, -Inf, Inf, -1e+7, 1e+7), nrow = 2)
ITCH <- function(par, dat, choice = FALSE, epsilon = 0.001, lambda = 0.001) {
  # Choice heuristic model function
  #
  # Args:
  #   par: Numeric vector of parameters.
  #   dat: Data frame with necessary columns G, R, D, and T.
  #   choice: Boolean indicating whether to return choice probability.
  #   epsilon: Numeric smoothing parameter between 0 and 1.
  #
  # Returns:
  #   If choice=TRUE, a numeric vector of choice probabilities.
  #   If choice=FALSE, a numeric representing the negative log-likelihood.
  
  # Compute the difference in option values
  diff <- (par[1] * dat$G + par[2] * dat$R + par[3] * dat$D + par[4] * dat$T)
  
  # Compute choice probability
  P <- 1 / (1 + exp(-par[5] * diff))
  
  # Apply smoothing
  P <- epsilon * 0.5 + (1 - epsilon) * P
  
  # If choice=TRUE, return probabilities
  if (choice) {
    return(P)
  } else {
    # Otherwise, compute and return negative log posterior
    P <- ifelse(as.logical(dat$LaterOptionChosen), P, 1 - P)
    nll <- -sum(log(P))
    nlpr <- -sum(dnorm(par, 0, 10, log = TRUE)) # negative log prior
    reg <- lambda * sum(par^2) # regularization term
    nlpo <- nll + nlpr + reg # negative log posterior
    return(nlpo)
  }
}

DRIFT_lim <- matrix(c(-Inf, Inf, -Inf, Inf, -Inf, Inf, -Inf, Inf, -1e+7, 1e+7), nrow = 2)
DRIFT <- function(par, dat, choice = FALSE, epsilon = 0.001, lambda = 0.001) {
  # DRIFT model function
  #
  # Args:
  #   par: Numeric vector of parameters.
  #   dat: Data frame with necessary columns DriftD, DriftR, DriftI, and DriftT.
  #   choice: Boolean indicating whether to return choice probability.
  #   epsilon: Numeric smoothing parameter between 0 and 1.
  #
  # Returns:
  #   If choice=TRUE, a numeric vector of choice probabilities.
  #   If choice=FALSE, a numeric representing the negative log-likelihood.
  
  # Compute the difference in option values
  diff <- (par[1] * dat$Drift_D + par[2] * dat$Drift_R + par[3] * dat$Drift_I + par[4] * dat$Drift_T)
  
  # Compute choice probability
  P <- 1 / (1 + exp(-par[5] * diff))
  
  # Apply smoothing
  P <- epsilon * 0.5 + (1 - epsilon) * P
  
  # If choice=TRUE, return probabilities
  if (choice) {
    return(P)
  } else {
    # Otherwise, compute and return negative log posterior
    P <- ifelse(as.logical(dat$LaterOptionChosen), P, 1 - P)
    nll <- -sum(log(P))
    nlpr <- -sum(dnorm(par, 0, 10, log = TRUE)) # negative log prior
    reg <- lambda * sum(par^2) # regularization term
    nlpo <- nll + nlpr + reg # negative log posterior
    return(nlpo)
  }
}

TRADE_lim <- matrix(c(1e-7, 1e+7, 1e-7, 1e+7, 1e-7, 1e+7, -1e+7, 1e+7), nrow = 2)
TRADE <- function(par, dat, choice = FALSE, epsilon = 0.001, lambda = 0.001) {
  # TRADE model function
  #
  # Args:
  #   par: Numeric vector of parameters.
  #   dat: Data frame with necessary columns X1, X2, T1, and T2.
  #   choice: Boolean indicating whether to return choice probability.
  #   epsilon: Numeric smoothing parameter between 0 and 1.
  #
  # Returns:
  #   If choice=TRUE, a numeric vector of choice probabilities.
  #   If choice=FALSE, a numeric representing the negative log-likelihood.
  
  # Utility conversion function
  .cnv <- function(x, g) {log(1 + g * x) / g}
  
  # Compute converted utilities
  a1 <- .cnv(dat$X2, par[2])
  a2 <- .cnv(dat$X1, par[2])
  a3 <- .cnv(dat$T2, par[3])
  a4 <- .cnv(dat$T1, par[3])
  
  # Compute the difference in option values
  diff <- ((a1 - a2) - par[1] * (a3 - a4))
  
  # Compute choice probability
  P <- 1 / (1 + exp(-par[4] * diff))
  
  # Apply smoothing
  P <- epsilon * 0.5 + (1 - epsilon) * P
  
  # If choice=TRUE, return probabilities
  if (choice) {
    return(P)
  } else {
    # Otherwise, compute and return negative log posterior
    P <- ifelse(as.logical(dat$LaterOptionChosen), P, 1 - P)
    nll <- -sum(log(P))
    nlpr <- -sum(dnorm(par, 0, 10, log = TRUE)) # negative log prior
    reg <- lambda * sum(par^2) # regularization term
    nlpo <- nll + nlpr + reg # negative log posterior
    return(nlpo)
  }
}
