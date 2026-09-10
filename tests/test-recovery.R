# Machinery tests: model sanity, keyed joins, and recovery on synthetic truth.
# Run from the project root with Rscript -e "testthat::test_dir('tests')".

library(testthat)
source(test_path("..", "R", "models.R"))
source(test_path("..", "R", "fit.R"))

cfg <- list(
  n_reps = 5L, n_restarts = 5L, train_prop = 0.7, epsilon = 1e-3, prior_sd = 10, seed = 1L
)
items <- readRDS(test_path("..", "data", "prepared", "choices.rds"))
items <- items[!duplicated(items[c("X1", "T1", "X2", "T2")]), ]

truth <- list(
  EXPO = c(delta = 0.7, sens = 3),
  HYPER2 = c(k = 0.5, s = 1, sens = 3),
  DEXPO = c(delta_1 = 0.3, delta_2 = 0.9, omega = 0.4, sens = 3),
  ITCH = c(beta_G = 1, beta_R = 0.5, beta_D = -1, beta_T = -0.5),
  DRIFT = c(beta_D = 1, beta_R = 0.5, beta_I = 0.5, beta_T = -1),
  TRADE = c(kappa = 1, gamma_x = 1, gamma_t = 1, sens = 50)
)

sample_trials <- function(subject, n, seed) {
  withr::with_seed(seed, {
    out <- items[sample.int(nrow(items), n, replace = TRUE), ]
    out$subject <- subject
    out$trial <- seq_len(n)
    out
  })
}

test_that("every model gives 0.5 for identical options and prefers a larger later amount", {
  same <- data.frame(
    X1 = 0.5, T1 = 1, X2 = 0.5, T2 = 1,
    itch_G = 0, itch_R = 0, itch_D = 0, itch_T = 0,
    drift_D = 0, drift_R = 0, drift_I = 0, drift_T = 0
  )
  for (model in names(MODELS)) {
    expect_equal(unname(MODELS[[model]]$prob(truth[[model]], same)), 0.5, info = model)
  }
  small <- items[1, ]
  large <- small
  large$X2 <- large$X2 * 2
  large$itch_G <- large$itch_G + 1
  large$drift_D <- large$drift_D + 1
  for (model in names(MODELS)) {
    p <- unname(MODELS[[model]]$prob(truth[[model]], rbind(small, large)))
    expect_true(all(p > 0 & p < 1), info = model)
    expect_gt(p[2], p[1], label = model)
  }
})

test_that("random starts stay inside the start box without warnings", {
  for (model in names(MODELS)) {
    spec <- MODELS[[model]]
    starts <- expect_no_warning(replicate(50, random_start(spec)))
    expect_true(all(starts >= spec$start_lower & starts <= spec$start_upper), info = model)
  }
})

test_that("the sample rule counts trials and minority choices per subject", {
  trials <- data.frame(
    subject = rep(c(7L, 3L), each = 5), later = c(1, 0, 0, 0, 0, 1, 1, 1, 0, 0)
  )
  flagged <- flag_subjects(trials)
  expect_equal(unique(flagged$n_trials), 5L)
  expect_equal(flagged$n_minority[flagged$subject == 7L], rep(1L, 5))
  expect_equal(flagged$n_minority[flagged$subject == 3L], rep(2L, 5))
  expect_equal(flagged$included, flagged$n_minority >= MIN_MINORITY)
})

test_that("splits are the right size and reproducible from the seed", {
  a <- withr::with_seed(3, make_splits(25, 4, 0.7))
  b <- withr::with_seed(3, make_splits(25, 4, 0.7))
  expect_identical(a, b)
  expect_true(all(lengths(a) == 17))
  expect_error(make_splits(2, 1, 0.7))
})

test_that("recovery metrics join by id, ignore row order and refuse unmatched rows", {
  true <- data.frame(
    subject = c(101L, 14L, 3L, 22L), model = "EXPO", parameter = "delta",
    value = c(0.2, 0.4, 0.6, 0.8), stringsAsFactors = FALSE
  )
  recovered <- true
  recovered$value <- true$value + 0.01
  shuffled <- recovered[c(3, 1, 4, 2), ]
  expect_equal(recovery_metrics(true, shuffled), recovery_metrics(true, recovered))
  expect_gt(recovery_metrics(true, shuffled)$pearson, 0.99)
  expect_lt(cor(true$value, shuffled$value), 0.9)
  expect_equal(recovery_metrics(true, recovered[-1, ])$n, 3)
  unmatched <- recovered
  unmatched$subject[1] <- 999L
  expect_error(recovery_metrics(true, unmatched))
})

test_that("the optimiser does at least as well as the truth on data from each model", {
  for (model in names(MODELS)) {
    trials <- sample_trials(1L, 200, seed = 11)
    trials <- simulate_trials(model, truth[[model]], trials, cfg, seed = 12)
    fit <- withr::with_seed(13, fit_model(model, trials, 10, cfg$epsilon, cfg$prior_sd))
    at_truth <- penalised_nll(truth[[model]], model, trials, cfg$epsilon, cfg$prior_sd)
    expect_lte(fit$value, at_truth + 1e-6, label = model)
    expect_named(fit$par, MODELS[[model]]$par_names)
  }
})

test_that("exponential agents are recovered from 400 trials", {
  n_agents <- 16
  agents <- do.call(rbind, lapply(seq_len(n_agents), function(i) {
    par <- withr::with_seed(20 + i, c(delta = runif(1, 0.3, 0.95), sens = runif(1, 1, 5)))
    data.frame(subject = i, model = "EXPO", parameter = names(par), value = unname(par))
  }))
  synthetic <- do.call(rbind, lapply(seq_len(n_agents), function(i) {
    simulate_trials(
      "EXPO", agent_par(agents, i, "EXPO"), sample_trials(i, 400, 30 + i), cfg, 40 + i
    )
  }))
  future::plan(future::sequential)
  recovered <- fit_agents(synthetic, "EXPO", cfg)
  metrics <- recovery_metrics(agents, recovered)
  expect_true(
    all(metrics$pearson > 0.9),
    info = paste(capture.output(print(metrics)), collapse = "\n")
  )
})

test_that("cross-validation picks the generating model for a separable pair", {
  synthetic <- do.call(rbind, lapply(1:4, function(i) {
    simulate_trials("EXPO", truth$EXPO, sample_trials(i, 100, 50 + i), cfg, 60 + i)
  }))
  future::plan(future::sequential)
  cv <- mccv_all(synthetic, c("EXPO", "TRADE"), cfg)
  expect_equal(nrow(cv), 4 * cfg$n_reps * 2)
  mean_loss <- aggregate(logloss ~ subject + model, cv, mean)
  winners <- vapply(
    split(mean_loss, mean_loss$subject), function(s) s$model[which.min(s$logloss)], ""
  )
  expect_gte(sum(winners == "EXPO"), 3)
})
