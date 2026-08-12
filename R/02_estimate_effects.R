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
  caliper = 0.20,
  std.caliper = TRUE,
  replace = FALSE
)

matched_data <- match_data(matched_design)
balance_object <- bal.tab(
  matched_design,
  un = TRUE,
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

clustered_program_effect <- function(model, data, label) {
  variance <- vcovCL(model, cluster = data$site_id, type = "HC1")
  estimate <- unname(coef(model)["program"])
  std_error <- unname(sqrt(diag(variance))["program"])

  tibble(
    estimator = label,
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error,
    n = nobs(model)
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
    "Matched ATT"
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
if (!is.finite(max_matched_smd) || max_matched_smd >= 0.10) {
  warning("At least one matched standardized difference is 0.10 or greater.")
}

print(effect_estimates)
