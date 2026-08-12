# Generate a reproducible, multisite administrative dataset.
# Treatment assignment intentionally depends on observed baseline factors so
# that a naive outcome comparison is confounded.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(20260812)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

n_sites <- 60L
n_people <- 2400L

site_frame <- tibble(
  site_id = sprintf("S%02d", seq_len(n_sites)),
  region = rep(c("North", "South", "West"), length.out = n_sites),
  site_outcome_effect = rnorm(n_sites, mean = 0, sd = 2.2)
)

evaluation_data <- tibble(
  participant_id = sprintf("P%04d", seq_len(n_people)),
  site_id = sample(site_frame$site_id, n_people, replace = TRUE),
  age = pmin(pmax(round(rnorm(n_people, 42, 12)), 18), 75),
  female = rbinom(n_people, 1, 0.56),
  urban = rbinom(n_people, 1, 0.62),
  prior_service = rbinom(n_people, 1, 0.31),
  baseline_score = pmin(pmax(rnorm(n_people, 50, 10), 15), 85)
) |>
  left_join(site_frame, by = "site_id") |>
  mutate(
    treatment_log_odds =
      -0.35 +
      0.65 * prior_service +
      0.35 * urban +
      0.18 * female -
      0.045 * (baseline_score - 50) +
      if_else(region == "South", 0.22, 0) -
      if_else(region == "West", 0.12, 0),
    treatment_probability = plogis(treatment_log_odds),
    program = rbinom(n(), 1, treatment_probability),
    outcome =
      18 +
      0.66 * baseline_score -
      0.05 * (age - 42) +
      0.90 * prior_service +
      2.75 * program +
      site_outcome_effect +
      rnorm(n(), mean = 0, sd = 6.5)
  ) |>
  transmute(
    participant_id,
    site_id,
    region,
    age,
    female,
    urban,
    prior_service,
    baseline_score = round(baseline_score, 2),
    program,
    outcome = round(outcome, 2)
  )

stopifnot(
  nrow(evaluation_data) == n_people,
  !anyDuplicated(evaluation_data$participant_id),
  all(evaluation_data$program %in% c(0, 1)),
  all(is.finite(evaluation_data$outcome))
)

write_csv(evaluation_data, "data/synthetic_program_data.csv")
message("Wrote ", nrow(evaluation_data), " synthetic records.")
