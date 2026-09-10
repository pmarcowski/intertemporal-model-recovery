# Computes every result the report reads.
#
# Stages: cross-validation of all models on the empirical subjects, agents (each
# model fit to each included subject), one recovery run per generating model, and
# a summary. Heavy stage output lands in output/raw/<stage>.rds tagged
# with a hash of the code, the configuration and the prepared data that produced
# it; a stage whose file carries the current hash is reused, so an interrupted run
# resumes and a code change recomputes. The compact summary the report reads is
# written to output/summary.rds. Run from the project root after the data
# preparation script.
#
# Environment overrides, for smoke runs: CMR_N_REPS (repetitions), CMR_N_SUBJECTS
# (first n included subjects; 0 means all), CMR_WORKERS (parallel workers).

suppressPackageStartupMessages({
  library(future)
  library(furrr)
})
source("R/models.R")
source("R/fit.R")

env_int <- function(name, default) {
  value <- Sys.getenv(name, "")
  if (value == "") return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 0) {
    stop(sprintf("%s must be a non-negative integer, got '%s'", name, value))
  }
  parsed
}

cfg <- list(
  n_reps = env_int("CMR_N_REPS", 100L),
  n_subjects = env_int("CMR_N_SUBJECTS", 0L),
  workers = env_int("CMR_WORKERS", max(1L, parallelly::availableCores() - 1L)),
  n_restarts = 10L,
  train_prop = 0.7,
  epsilon = 1e-3,
  prior_sd = 10,
  seed = 42L
)
stopifnot(cfg$n_reps >= 1, cfg$workers >= 1)

trials <- readRDS("data/prepared/choices.rds")
included <- trials[trials$included, ]
if (cfg$n_subjects > 0) {
  keep <- head(unique(included$subject), cfg$n_subjects)
  trials <- trials[trials$subject %in% keep, ]
  included <- included[included$subject %in% keep, ]
}
models <- names(MODELS)

run_hash <- digest::digest(
  list(
    readLines("R/models.R"), readLines("R/fit.R"),
    cfg[setdiff(names(cfg), "workers")], digest::digest(trials)
  ),
  algo = "sha256"
)

save_rds <- function(x, path) {
  # Writes through a temporary file and retries, because a file under sync can be
  # locked for a few seconds.
  tmp <- paste0(path, ".tmp")
  last_error <- "rename returned FALSE"
  for (attempt in 1:3) {
    ok <- tryCatch(
      {
        saveRDS(x, tmp)
        file.rename(tmp, path)
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        FALSE
      }
    )
    if (isTRUE(ok)) return(invisible(path))
    if (attempt < 3) Sys.sleep(10)
  }
  stop(sprintf("could not write %s in 3 attempts; last error: %s", path, last_error))
}

stage <- function(name, compute) {
  path <- file.path("output/raw", paste0(name, ".rds"))
  if (file.exists(path)) {
    previous <- readRDS(path)
    if (identical(attr(previous, "run_hash"), run_hash)) {
      message(sprintf("%-16s reused", name))
      return(previous)
    }
  }
  started <- Sys.time()
  result <- compute()
  attr(result, "run_hash") <- run_hash
  attr(result, "elapsed_s") <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  save_rds(result, path)
  message(sprintf("%-16s %7.1f s", name, attr(result, "elapsed_s")))
  result
}

# Summary helpers ---------------------------------------------------------

mean_logloss <- function(cv) {
  # Per-subject mean held-out log-loss of each model over repetitions.
  out <- aggregate(logloss ~ subject + model, cv, mean)
  out <- out[order(out$subject, out$model), ]
  rownames(out) <- NULL
  out
}

subject_table <- function(trials) {
  columns <- c("subject", "n_trials", "n_later", "n_minority", "included")
  out <- trials[!duplicated(trials$subject), columns]
  rownames(out) <- NULL
  out
}

determinism <- function(generating, agents, trials) {
  # Mean |P(later) - 0.5| of each agent over its own trials: how far from coin
  # flips the generating model's choices are.
  by_subject <- split(trials, trials$subject)
  vapply(by_subject, function(dat) {
    par <- agent_par(agents, dat$subject[1], generating)
    mean(abs(choice_prob(generating, par, dat, cfg$epsilon) - 0.5))
  }, numeric(1))
}

summarise_run <- function(trials, included, empirical_cv, agents, recovery) {
  synthetic_subjects <- do.call(rbind, lapply(recovery, function(r) {
    out <- subject_table(r$synthetic)
    out$determinism <- determinism(r$generating, agents, included)[as.character(out$subject)]
    cbind(generating = r$generating, out, stringsAsFactors = FALSE)
  }))
  rownames(synthetic_subjects) <- NULL
  recovery_cv <- do.call(rbind, lapply(recovery, function(r) {
    cbind(generating = r$generating, mean_logloss(r$cv), stringsAsFactors = FALSE)
  }))
  rownames(recovery_cv) <- NULL
  recovered <- do.call(rbind, lapply(recovery, function(r) r$recovered))
  rownames(recovered) <- NULL
  list(
    config = cfg[setdiff(names(cfg), "workers")],
    run_hash = run_hash,
    subjects = subject_table(trials),
    empirical = mean_logloss(empirical_cv),
    agents = agents,
    synthetic_subjects = synthetic_subjects,
    recovery = recovery_cv,
    recovered = recovered
  )
}

# Run ---------------------------------------------------------------------

message(sprintf(
  "%d subjects (%d included), %d models, %d repetitions, %d workers",
  length(unique(trials$subject)), length(unique(included$subject)),
  length(models), cfg$n_reps, cfg$workers
))
dir.create("output/raw", recursive = TRUE, showWarnings = FALSE)
plan(multisession, workers = cfg$workers)

empirical_cv <- stage("empirical_cv", function() mccv_all(trials, models, cfg))
agents <- stage("agents", function() fit_agents(included, models, cfg))
recovery <- lapply(models, function(generating) {
  stage(
    paste0("recovery_", generating),
    function() recover_generating(generating, agents, included, models, cfg)
  )
})
names(recovery) <- models

plan(sequential)
summary <- summarise_run(trials, included, empirical_cv, agents, recovery)
save_rds(summary, "output/summary.rds")
message("summary written to output/summary.rds")
