# Estimate the ATT with matching and two prespecified robustness models.

suppressPackageStartupMessages({
  library(cobalt)
  library(dplyr)
  library(lmtest)
  library(MatchIt)
  library(readr)
  library(sandwich)
  library(tibble)
})

input_path <- "data/synthetic_program_data.csv"
if (!file.exists(input_path)) {
  stop("Run R/01_generate_synthetic_data.R before estimating effects.")
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
evaluation_data <- read_csv(input_path, show_col_types = FALSE)

propensity_formula <- program ~
  baseline_score + age + female + urban + prior_service + factor(region)

matched_design <- matchit(
  formula = propensity_formula,
  data = evaluation_data,
  method = "nearest",
  distance = "glm",
  estimand = "ATT",
  ratio = 1,
  # Locked after an outcome-blind balance audit; see docs/design-revision.md.
  caliper = 0.10,
  std.caliper = TRUE,
  replace = FALSE
)

matched_data <- match_data(matched_design)
balance_object <- bal.tab(
  matched_design,
  un = TRUE,
  binary = "std",
  thresholds = c(m = 0.10)
)

balance_table <- balance_object$Balance |>
  as.data.frame() |>
  rownames_to_column("covariate") |>
  select(any_of(c("covariate", "Type", "Diff.Un", "Diff.Adj"))) |>
  rename(
    variable_type = Type,
    smd_unadjusted = Diff.Un,
    smd_matched = Diff.Adj
  )

write_csv(balance_table, "outputs/balance_table.csv")

matching_flow <- bind_rows(lapply(0:1, function(group) {
  before <- sum(evaluation_data$program == group)
  after <- sum(matched_data$program == group)
  tibble(group = c("Comparison", "Program")[group + 1],
         before = before, matched = after, unmatched = before - after,
         retained_percent = 100 * after / before)
}))
write_csv(matching_flow, "outputs/matching_flow.csv")
matched_label <- if (any(matching_flow$unmatched[matching_flow$group == "Program"] > 0)) {
  "Matched-treated subset"
} else "Matched ATT"

clustered_program_effect <- function(model, data, label) {
  variance <- vcovCL(model, cluster = data$site_id, type = "HC1")
  estimate <- unname(coef(model)["program"])
  std_error <- unname(sqrt(diag(variance))["program"])
  clusters <- length(unique(data$site_id))
  critical <- qt(.975, df = clusters - 1)

  tibble(
    estimator = label,
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - critical * std_error,
    conf_high = estimate + critical * std_error,
    n = nobs(model),
    clusters = clusters,
    degrees_freedom = clusters - 1
  )
}

unadjusted_model <- lm(outcome ~ program, data = evaluation_data)
adjusted_model <- lm(
  outcome ~ program + baseline_score + age + female + urban +
    prior_service + factor(region),
  data = evaluation_data
)
matched_model <- lm(
  outcome ~ program + baseline_score + age + female + urban +
    prior_service + factor(region),
  data = matched_data,
  weights = weights
)

propensity_model <- glm(
  propensity_formula,
  data = evaluation_data,
  family = binomial()
)
evaluation_data <- evaluation_data |>
  mutate(
    propensity_score = pmin(
      pmax(predict(propensity_model, type = "response"), 0.01),
      0.99
    ),
    att_weight = if_else(
      program == 1,
      1,
      propensity_score / (1 - propensity_score)
    )
  )
weighted_model <- lm(
  outcome ~ program + baseline_score + age + female + urban +
    prior_service + factor(region),
  data = evaluation_data,
  weights = att_weight
)

effect_estimates <- bind_rows(
  clustered_program_effect(
    unadjusted_model,
    evaluation_data,
    "Unadjusted difference"
  ),
  clustered_program_effect(
    adjusted_model,
    evaluation_data,
    "Regression adjusted"
  ),
  clustered_program_effect(
    matched_model,
    matched_data,
    matched_label
  ),
  clustered_program_effect(
    weighted_model,
    evaluation_data,
    "ATT weighted"
  )
) |>
  mutate(
    across(
      c(estimate, std_error, conf_low, conf_high),
      ~ round(.x, 3)
    )
  )

write_csv(effect_estimates, "outputs/effect_estimates.csv")

max_matched_smd <- max(abs(balance_table$smd_matched), na.rm = TRUE)
stopifnot(is.finite(max_matched_smd), max_matched_smd < 0.10)

write_csv(tibble(group = c("Program", "Comparison"),
  effective_n = sapply(1:0, function(g) {
    w <- evaluation_data$att_weight[evaluation_data$program == g]
    sum(w)^2 / sum(w^2)
  }), maximum_weight = sapply(1:0, function(g) max(evaluation_data$att_weight[evaluation_data$program == g]))),
  "outputs/weight_diagnostics.csv")

print(effect_estimates)
