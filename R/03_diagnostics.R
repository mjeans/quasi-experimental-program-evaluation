# Create design diagnostics after the estimators have been run.

suppressPackageStartupMessages({
  library(cobalt)
  library(dplyr)
  library(ggplot2)
  library(MatchIt)
  library(readr)
})

input_path <- "data/synthetic_program_data.csv"
source("R/portfolio_report.R")
if (!file.exists(input_path)) {
  stop("Run R/01_generate_synthetic_data.R before creating diagnostics.")
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
evaluation_data <- read_csv(input_path, show_col_types = FALSE)

propensity_formula <- program ~
  baseline_score + age + female + urban + prior_service + factor(region)
propensity_model <- glm(
  propensity_formula,
  data = evaluation_data,
  family = binomial()
)
evaluation_data <- evaluation_data |>
  mutate(propensity_score = predict(propensity_model, type = "response"))

overlap_plot <- ggplot(
  evaluation_data,
  aes(
    x = propensity_score,
    fill = factor(program),
    color = factor(program)
  )
) +
  geom_density(alpha = 0.22, linewidth = 0.8) +
  scale_fill_manual(
    values = c("#64748B", "#15803D"),
    labels = c("Comparison", "Program")
  ) +
  scale_color_manual(
    values = c("#475569", "#166534"),
    labels = c("Comparison", "Program")
  ) +
  labs(
    title = "Propensity-score overlap before matching",
    x = "Estimated propensity score",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(
  "outputs/propensity_overlap.png",
  overlap_plot,
  width = 8,
  height = 5,
  dpi = 160
)

matched_design <- matchit(
  propensity_formula,
  data = evaluation_data,
  method = "nearest",
  distance = "glm",
  estimand = "ATT",
  ratio = 1,
  caliper = 0.10,
  std.caliper = TRUE
)

png("outputs/love_plot.png", width = 1400, height = 900, res = 160)
love.plot(
  matched_design,
  stats = "mean.diffs",
  binary = "std",
  abs = TRUE,
  thresholds = c(m = 0.10),
  var.order = "unadjusted",
  colors = c("#64748B", "#15803D"),
  shapes = c(16, 17),
  sample.names = c("Before matching", "After matching")
)
dev.off()

balance_plot <- love.plot(matched_design, stats = "mean.diffs", binary = "std",
  abs = TRUE, thresholds = c(m = .10), var.order = "unadjusted",
  colors = c("#687b89", "#087e83"), shapes = c(16, 17),
  sample.names = c("Before matching", "After matching")) +
  labs(title = "Covariate balance before and after matching",
       subtitle = "Synthetic multisite evaluation; all differences standardized")
publish_plot(overlap_plot, "assets/propensity-overlap.svg", "Propensity-score overlap",
             "Synthetic propensity-score distributions before matching; overlap is necessary but not sufficient for causal identification.")
publish_plot(balance_plot, "assets/covariate-balance.svg", "Covariate balance",
             "Absolute standardized differences before and after matching; reference line at 0.10.")
effects <- read_csv("outputs/effect_estimates.csv", show_col_types = FALSE)
effects$estimator <- factor(effects$estimator, levels = rev(effects$estimator))
effect_plot <- ggplot(effects, aes(estimate, estimator)) +
  geom_vline(xintercept = 0, color = "#687b89", linetype = 2) +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = estimator), color = "#087e83", linewidth = 1) +
  geom_point(color = "#087e83", size = 3) +
  labs(title = "Effect estimates depend on the design",
       subtitle = "Synthetic data; site-clustered HC1 variance and t intervals (G - 1 df)",
       x = "Outcome-point difference (95% CI)", y = NULL,
       caption = "Matched-treated subset and full treated population are not interchangeable targets.")
publish_plot(effect_plot, "assets/effect-estimates.svg", "Evaluation estimates with uncertainty",
             "Four descriptive or adjusted estimates; intervals use site-clustered standard errors and a t reference.")
flow <- read_csv("outputs/matching_flow.csv", show_col_types = FALSE)
write_research_report(c("# Executed program-evaluation report", "",
  "## Question and design", "",
  "Can a prespecified matching/weighting workflow recover the constant 2.75-point effect in a synthetic multisite evaluation? Treatment depends on observed baseline factors. Real-world identification additionally requires no unmeasured confounding, positivity, and consistency.", "",
  "## Who remains matched?", "", md_table(flow), "",
  "The standardized propensity-score caliper is 0.10. It was tightened after an outcome-blind balance audit; see [design revision](../docs/design-revision.md). The original 0.20 rule failed the declared 0.10 balance threshold.", "",
  "A caliper can exclude treated observations. The matched analysis targets the retained treated subset when this occurs; weighting targets the full treated population under its assumptions.", "",
  "## Design diagnostics", "", "![Covariate balance](../assets/covariate-balance.svg)", "",
  "![Propensity overlap](../assets/propensity-overlap.svg)", "",
  "[Weight diagnostics](weight_diagnostics.csv). Propensity scores are clipped to 0.01–0.99 before ATT weighting; this numerical rule is not a substitute for overlap assessment.", "",
  "## Estimates", "", "![Effect estimates](../assets/effect-estimates.svg)", "", md_table(effects), "",
  "## Limitations", "", "These are synthetic recovery checks, not evidence of a real program effect. Site-clustered HC1 intervals do not capture every source of matching/propensity-estimation uncertainty. Caliper restrictions change the analyzed population. No unmeasured-confounding sensitivity model is claimed."))
