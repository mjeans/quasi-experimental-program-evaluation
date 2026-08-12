# Quasi-experimental program evaluation

A reproducible evaluation of a hypothetical service program using synthetic administrative data. The project estimates the average treatment effect on participants while making the design assumptions, balance diagnostics, and sensitivity of the estimate visible.

This is a portfolio demonstration. The data are generated in the repository and do not represent a real organization, client, or participant.

## Evaluation question

Among eligible participants, what was the effect of program participation on the follow-up outcome after accounting for observed baseline differences between participants and nonparticipants?

## Design

1. Generate a synthetic multisite administrative dataset with realistic selection into the program.
2. Estimate propensity scores from prespecified baseline covariates.
3. Perform 1:1 nearest-neighbor matching with a caliper.
4. Assess standardized mean differences before and after matching.
5. Estimate the participant effect with site-clustered standard errors.
6. Compare the matched estimate with adjusted-regression and ATT-weighted estimates.

The estimand is the average treatment effect on the treated (ATT). The workflow treats the propensity model as a design tool: outcome results are not interpreted until common support and balance have been reviewed.

## Repository structure

~~~text
R/
  01_generate_synthetic_data.R
  02_estimate_effects.R
  03_diagnostics.R
analysis/
  program_evaluation.qmd
stata/
  quasi_experimental_analysis.do
data/
  README.md
.github/workflows/
  validate-analysis.yml
~~~

## Run the analysis

With R 4.4 or later:

~~~r
install.packages(c(
  "cobalt", "dplyr", "ggplot2", "lmtest",
  "MatchIt", "readr", "sandwich", "tibble"
))

source("R/01_generate_synthetic_data.R")
source("R/02_estimate_effects.R")
source("R/03_diagnostics.R")
~~~

The companion Stata workflow reads the same generated CSV and estimates matching and doubly robust models. It is included to demonstrate equivalent implementation across the two environments; GitHub Actions validates the R workflow because Stata requires a commercial license.

## Outputs

Running the project creates:

- `outputs/effect_estimates.csv` with unadjusted, regression-adjusted, matched, and ATT-weighted estimates
- `outputs/balance_table.csv` with pre/post-match standardized mean differences
- `outputs/love_plot.png` for covariate balance
- `outputs/propensity_overlap.png` for common-support assessment

## Interpretation boundaries

Matching can improve comparability on observed baseline variables, but it cannot remove bias from unmeasured confounding. A real evaluation would also document intervention timing, attrition, data provenance, spillovers, missing-data mechanisms, and the plausibility of conditional exchangeability with subject-matter experts.

## Skills demonstrated

Program evaluation · causal inference · propensity-score matching · balance diagnostics · clustered inference · reproducible analysis · R · Stata

## License

MIT
