effects <- read.csv("outputs/effect_estimates.csv")
balance <- read.csv("outputs/balance_table.csv")
flow <- read.csv("outputs/matching_flow.csv")
stopifnot(nrow(effects) == 4L, all(is.finite(effects$estimate)),
          all(effects$conf_low < effects$estimate), all(effects$conf_high > effects$estimate),
          all(effects$clusters == 60), all(effects$degrees_freedom == 59),
          max(abs(balance$smd_matched)) < .10,
          all(flow$before == flow$matched + flow$unmatched))
# A prespecified tolerance for this fixed synthetic data-generating process.
stopifnot(all(abs(effects$estimate[effects$estimator != "Unadjusted difference"] - 2.75) < 1.0))
cat("Balance, matching flow, clustered intervals, and synthetic recovery checks passed.\n")
