# Choice models for binary intertemporal choice.
#
# Every trial offers a smaller-sooner amount X1 at delay T1 against a larger-later
# amount X2 at delay T2. Amounts arrive scaled to the largest later amount and the
# attribute columns used by the heuristic models arrive standardised, both from
# R/prepare_data.R. Each entry of MODELS supplies the probability of choosing the
# later option before lapse smoothing, the box constraints used in fitting, the
# box (and scale) its random starts are drawn from, and parameter names.
#
# Discounting models value each option and choose by a power rule; heuristic
# models weigh trial attributes and choose by a logistic rule whose weights carry
# the scale, so they have no separate sensitivity parameter. Definitions follow
# Wulff and van den Bos (2018).
#
# The choice rules and the registry share one local environment so that every
# probability function carries its helpers when it is shipped to a parallel
# worker.

MODELS <- local({
  power_rule <- function(u_later, u_sooner, sens) {
    u_later^sens / (u_later^sens + u_sooner^sens)
  }

  logistic_rule <- function(diff) {
    1 / (1 + exp(-diff))
  }

  list(
    EXPO = list(
      label = "Exponential",
      family = "discounting",
      par_names = c("delta", "sens"),
      lower = c(1e-3, -10),
      upper = c(100, 10),
      start_lower = c(1e-3, -10),
      start_upper = c(100, 10),
      log_start = c(TRUE, FALSE),
      prob = function(par, dat) {
        power_rule(dat$X2 * par[1]^dat$T2, dat$X1 * par[1]^dat$T1, par[2])
      }
    ),
    HYPER2 = list(
      label = "Hyperboloid",
      family = "discounting",
      par_names = c("k", "s", "sens"),
      lower = c(1e-3, 1e-7, -10),
      upper = c(10, 10, 10),
      start_lower = c(1e-3, 1e-7, -10),
      start_upper = c(10, 10, 10),
      log_start = c(TRUE, TRUE, FALSE),
      prob = function(par, dat) {
        power_rule(
          dat$X2 / (1 + par[1] * dat$T2)^par[2],
          dat$X1 / (1 + par[1] * dat$T1)^par[2],
          par[3]
        )
      }
    ),
    DEXPO = list(
      label = "Dual exponential",
      family = "discounting",
      par_names = c("delta_1", "delta_2", "omega", "sens"),
      lower = c(1e-3, 1e-3, 1e-7, -10),
      upper = c(1, 1, 1, 10),
      start_lower = c(1e-3, 1e-3, 1e-7, -10),
      start_upper = c(1, 1, 1, 10),
      log_start = c(TRUE, TRUE, TRUE, FALSE),
      prob = function(par, dat) {
        weight <- function(t) par[3] * par[1]^t + (1 - par[3]) * par[2]^t
        power_rule(dat$X2 * weight(dat$T2), dat$X1 * weight(dat$T1), par[4])
      }
    ),
    ITCH = list(
      label = "Intertemporal choice heuristic",
      family = "heuristic",
      par_names = c("beta_G", "beta_R", "beta_D", "beta_T"),
      lower = rep(-Inf, 4),
      upper = rep(Inf, 4),
      start_lower = rep(-3, 4),
      start_upper = rep(3, 4),
      log_start = rep(FALSE, 4),
      prob = function(par, dat) {
        logistic_rule(
          par[1] * dat$itch_G + par[2] * dat$itch_R + par[3] * dat$itch_D + par[4] * dat$itch_T
        )
      }
    ),
    DRIFT = list(
      label = "DRIFT",
      family = "heuristic",
      par_names = c("beta_D", "beta_R", "beta_I", "beta_T"),
      lower = rep(-Inf, 4),
      upper = rep(Inf, 4),
      start_lower = rep(-3, 4),
      start_upper = rep(3, 4),
      log_start = rep(FALSE, 4),
      prob = function(par, dat) {
        logistic_rule(
          par[1] * dat$drift_D + par[2] * dat$drift_R + par[3] * dat$drift_I + par[4] * dat$drift_T
        )
      }
    ),
    TRADE = list(
      label = "Trade-off",
      family = "heuristic",
      par_names = c("kappa", "gamma_x", "gamma_t", "sens"),
      lower = c(1e-7, 1e-7, 1e-7, -1e7),
      upper = c(1e7, 1e7, 1e7, 1e7),
      start_lower = c(1e-3, 1e-3, 1e-3, 1),
      start_upper = c(1e3, 1e3, 1e3, 1e5),
      log_start = c(TRUE, TRUE, TRUE, TRUE),
      prob = function(par, dat) {
        convert <- function(x, gamma) log1p(gamma * x) / gamma
        gain <- convert(dat$X2, par[2]) - convert(dat$X1, par[2])
        wait <- convert(dat$T2, par[3]) - convert(dat$T1, par[3])
        logistic_rule(par[4] * (gain - par[1] * wait))
      }
    )
  )
})

choice_prob <- function(model, par, dat, epsilon) {
  # Probability of choosing the later option with a symmetric lapse of size epsilon.
  epsilon / 2 + (1 - epsilon) * MODELS[[model]]$prob(par, dat)
}

penalised_nll <- function(par, model, dat, epsilon, prior_sd) {
  # Negative log-likelihood of the observed choices plus a weak Gaussian prior on
  # every parameter, the objective minimised in fitting.
  p <- choice_prob(model, par, dat, epsilon)
  -sum(ifelse(dat$later == 1L, log(p), log1p(-p))) - sum(dnorm(par, 0, prior_sd, log = TRUE))
}
