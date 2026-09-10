# Prepares the analysis sample from the raw Ericson et al. (2015) export.
#
# Reads data/raw/choices.csv and writes data/prepared/choices.rds with one row
# per answered trial of condition 1: integer subject and trial ids, the two
# options, the choice, the standardised attributes the heuristic models use, and
# the per-subject counts behind the inclusion flag. Amounts are scaled to the
# largest later amount so the discounting models see values in (0, 1].
# Run from the project root.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("R/models.R")
source("R/fit.R")

raw <- read_csv("data/raw/choices.csv", col_types = cols(.default = col_double()))
needed <- c("Subject", "Condition", "Question", "X1", "T1", "X2", "T2", "LaterOptionChosen")
stopifnot(all(needed %in% names(raw)))

trials <- raw %>%
  filter(Condition == 1) %>%
  transmute(
    subject = as.integer(Subject), trial = as.integer(Question),
    X1, T1, X2, T2, later = as.integer(LaterOptionChosen)
  ) %>%
  filter(!is.na(later)) %>%
  group_by(subject) %>%
  filter(n() >= MIN_TRIALS) %>%
  ungroup() %>%
  arrange(subject, trial)

stopifnot(
  all(trials$X1 < trials$X2), all(trials$T1 < trials$T2),
  !anyDuplicated(trials[c("subject", "trial")]), all(trials$later %in% 0:1)
)

standardise <- function(x) as.numeric(scale(x))
x_max <- max(trials$X2)

trials <- trials %>%
  mutate(
    x_mid = (X1 + X2) / 2, t_mid = (T1 + T2) / 2,
    itch_G = standardise(X2 - x_mid),
    itch_R = standardise((X2 - X1) / x_mid),
    itch_D = standardise(T2 - t_mid),
    itch_T = standardise((T2 - T1) / t_mid),
    drift_D = standardise(X2 - X1),
    drift_R = standardise((X2 - X1) / X1),
    drift_I = standardise((X2 / X1)^(1 / (T2 - T1)) - 1),
    drift_T = standardise(T2 - T1),
    X1 = X1 / x_max, X2 = X2 / x_max
  ) %>%
  select(-x_mid, -t_mid) %>%
  as.data.frame() %>%
  flag_subjects()

stopifnot(!anyNA(trials))
saveRDS(trials, "data/prepared/choices.rds")

subjects <- trials[!duplicated(trials$subject), ]
message(sprintf(
  paste(
    "condition 1: %d subjects in the raw export, %d with at least %d answered trials,",
    "%d included (at least %d choices of the rarer option)"
  ),
  n_distinct(raw$Subject[raw$Condition == 1]), nrow(subjects), MIN_TRIALS,
  sum(subjects$included), MIN_MINORITY
))
message(sprintf("%d trials written to data/prepared/choices.rds", nrow(trials)))
