# Sample rules, fitting, cross-validation, simulation and recovery machinery.
# Requires R/models.R. Everything stochastic derives from a seed passed in, so a
# result is identical whatever the number of workers.

MIN_TRIALS <- 20L   # answered trials a subject needs to enter the sample
MIN_MINORITY <- 3L  # choices of the rarer option a subject needs to be analysed

flag_subjects <- function(trials) {
  # Adds per-subject trial counts and the inclusion flag to a trials data frame.
  n_trials <- tapply(trials$later, trials$subject, length)
  n_later <- tapply(trials$later, trials$subject, sum)
  idx <- match(trials$subject, as.integer(names(n_trials)))
  trials$n_trials <- as.integer(n_trials[idx])
  trials$n_later <- as.integer(n_later[idx])
  trials$n_minority <- pmin(trials$n_later, trials$n_trials - trials$n_later)
  trials$included <- trials$n_minority >= MIN_MINORITY
  trials
}

# Fitting -----------------------------------------------------------------

random_start <- function(spec) {
  # One start vector drawn from the model's start box, log-uniform where flagged.
  lo <- spec$start_lower
  hi <- spec$start_upper
  u <- runif(length(lo))
  start <- lo + u * (hi - lo)
  lg <- spec$log_start
  start[lg] <- 10^(log10(lo[lg]) + u[lg] * (log10(hi[lg]) - log10(lo[lg])))
  start
}

fit_model <- function(model, trials, n_restarts, epsilon, prior_sd) {
  # Best of n_restarts L-BFGS-B fits from random starts. A start whose objective
  # is not finite is discarded and redrawn, up to three times n_restarts attempts.
  spec <- MODELS[[model]]
  best <- NULL
  n_ok <- 0L
  n_tried <- 0L
  last_error <- "objective not finite"
  while (n_ok < n_restarts && n_tried < 3L * n_restarts) {
    n_tried <- n_tried + 1L
    fit <- tryCatch(
      optim(
        random_start(spec), penalised_nll,
        model = model, dat = trials, epsilon = epsilon, prior_sd = prior_sd,
        method = "L-BFGS-B", lower = spec$lower, upper = spec$upper
      ),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(fit) || !is.finite(fit$value)) next
    n_ok <- n_ok + 1L
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  if (is.null(best)) {
    stop(sprintf(
      "no successful fit in %d attempts: model %s, subject %s, %d trials; last error: %s",
      n_tried, model, trials$subject[1], nrow(trials), last_error
    ))
  }
  list(par = setNames(best$par, spec$par_names), value = best$value, n_ok = n_ok)
}

test_logloss <- function(model, par, trials, epsilon) {
  # Mean negative log-probability of the held-out choices.
  p <- choice_prob(model, par, trials, epsilon)
  -mean(ifelse(trials$later == 1L, log(p), log1p(-p)))
}

# Cross-validation --------------------------------------------------------

make_splits <- function(n, n_reps, prop) {
  # Training-row indices for n_reps random splits of n trials.
  n_train <- floor(n * prop)
  stopifnot(n_train >= 2, n - n_train >= 1)
  replicate(n_reps, sort(sample.int(n, n_train)), simplify = FALSE)
}

mccv_subject <- function(trials, models, cfg, seed) {
  # Monte Carlo cross-validation for one subject. Every model is fit to the same
  # training split of each repetition and scored on the same held-out trials.
  withr::with_seed(seed, {
    splits <- make_splits(nrow(trials), cfg$n_reps, cfg$train_prop)
    rows <- vector("list", cfg$n_reps * length(models))
    k <- 0L
    for (r in seq_along(splits)) {
      train <- trials[splits[[r]], ]
      test <- trials[-splits[[r]], ]
      for (model in models) {
        fit <- fit_model(model, train, cfg$n_restarts, cfg$epsilon, cfg$prior_sd)
        k <- k + 1L
        rows[[k]] <- data.frame(
          subject = trials$subject[1], rep = r, model = model,
          n_train = nrow(train), n_test = nrow(test),
          logloss = test_logloss(model, fit$par, test, cfg$epsilon),
          n_ok = fit$n_ok,
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  })
}

mccv_all <- function(trials, models, cfg) {
  # Cross-validates every subject in trials, in parallel over subjects.
  by_subject <- split(trials, trials$subject)
  ids <- as.integer(names(by_subject))
  rows <- furrr::future_map2(
    by_subject, ids,
    function(dat, id) mccv_subject(dat, models, cfg, cfg$seed + id),
    .options = furrr::furrr_options(seed = TRUE)
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  stopifnot(nrow(out) == length(ids) * cfg$n_reps * length(models))
  out
}

# Agents and simulation ---------------------------------------------------

fit_agents <- function(trials, models, cfg) {
  # Fits each model to each subject's full data. Returns one row per subject,
  # model and parameter.
  by_subject <- split(trials, trials$subject)
  ids <- as.integer(names(by_subject))
  grid <- expand.grid(subject = ids, model = models, stringsAsFactors = FALSE)
  fits <- furrr::future_map2(
    grid$subject, grid$model,
    function(id, model) {
      fit <- withr::with_seed(
        cfg$seed + id,
        fit_model(model, by_subject[[as.character(id)]], cfg$n_restarts, cfg$epsilon, cfg$prior_sd)
      )
      data.frame(
        subject = id, model = model, parameter = names(fit$par), value = unname(fit$par),
        objective = fit$value, stringsAsFactors = FALSE
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  out <- do.call(rbind, fits)
  rownames(out) <- NULL
  n_par <- sum(vapply(models, function(m) length(MODELS[[m]]$par_names), 1L))
  stopifnot(nrow(out) == n_par * length(ids))
  out
}

agent_par <- function(agents, id, model) {
  # The fitted parameter vector of one agent, in the model's parameter order.
  rows <- agents[agents$subject == id & agents$model == model, ]
  stopifnot(nrow(rows) == length(MODELS[[model]]$par_names))
  setNames(rows$value, rows$parameter)[MODELS[[model]]$par_names]
}

simulate_trials <- function(model, par, trials, cfg, seed) {
  # Replaces the choice column of a subject's trials with draws from the model.
  withr::with_seed(seed, {
    trials$later <- rbinom(nrow(trials), 1L, choice_prob(model, par, trials, cfg$epsilon))
    trials
  })
}

# Recovery ----------------------------------------------------------------

recover_generating <- function(generating, agents, trials, models, cfg) {
  # Model and parameter recovery for one generating model: each agent simulates
  # one synthetic dataset on its own trials, the sample rule is reapplied, every
  # model is cross-validated on the synthetic subjects, and the generating model
  # is refit to them for parameter recovery.
  by_subject <- split(trials, trials$subject)
  ids <- as.integer(names(by_subject))
  offset <- 100000L * match(generating, names(MODELS))
  synthetic <- do.call(rbind, lapply(ids, function(id) {
    simulate_trials(
      generating, agent_par(agents, id, generating), by_subject[[as.character(id)]],
      cfg, cfg$seed + offset + id
    )
  }))
  rownames(synthetic) <- NULL
  synthetic <- flag_subjects(synthetic)
  analysed <- synthetic[synthetic$included, ]
  stopifnot(nrow(analysed) > 0)
  list(
    generating = generating,
    synthetic = synthetic,
    cv = mccv_all(analysed, models, cfg),
    recovered = fit_agents(analysed, generating, cfg)
  )
}

recovery_metrics <- function(true, recovered) {
  # Agreement between generating and recovered parameters, joined by subject,
  # model and parameter, never by position.
  joined <- merge(
    true, recovered, by = c("subject", "model", "parameter"),
    suffixes = c("_true", "_recovered")
  )
  stopifnot(nrow(joined) == nrow(recovered))
  groups <- split(joined, list(joined$model, joined$parameter), drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    data.frame(
      model = g$model[1], parameter = g$parameter[1], n = nrow(g),
      pearson = cor(g$value_true, g$value_recovered),
      spearman = cor(g$value_true, g$value_recovered, method = "spearman"),
      rmse = sqrt(mean((g$value_true - g$value_recovered)^2)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}
