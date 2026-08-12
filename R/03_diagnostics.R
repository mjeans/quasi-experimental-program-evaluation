# Create design diagnostics after the estimators have been run.

suppressPackageStartupMessages({
  library(cobalt)
  library(dplyr)
  library(ggplot2)
  library(MatchIt)
  library(readr)
})

input_path <- "data/synthetic_program_data.csv"
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
  caliper = 0.20,
  std.caliper = TRUE
)

png("outputs/love_plot.png", width = 1400, height = 900, res = 160)
love.plot(
  matched_design,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.10),
  var.order = "unadjusted",
  colors = c("#64748B", "#15803D"),
  shapes = c(16, 17),
  sample.names = c("Before matching", "After matching")
)
dev.off()
