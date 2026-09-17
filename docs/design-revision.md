# Outcome-blind matching design revision

The September 2026 audit found that the original standardized caliper of 0.20 did not satisfy the repository's declared absolute standardized mean difference (SMD) threshold of 0.10. Binary covariates are now standardized consistently, too.

Using only treatment and baseline covariates, the fixed synthetic sample was checked at three calipers with the same nearest-neighbor algorithm, no replacement, and 1:1 matching:

| Standardized caliper | Largest absolute post-match SMD | Treated retained (of 1,301) |
|---|---:|---:|
| 0.20 (original) | 0.1594 | 974 |
| 0.10 (selected) | 0.0622 | 918 |
| 0.05 | 0.0398 | 898 |

The least restrictive checked rule meeting the balance threshold, 0.10, is now fixed in both R scripts. Outcomes and estimated effects were not used to choose it. This is a documented diagnostic-driven revision, not a claim that 0.10 was prespecified before this audit. A future independent study should prespecify its design-review procedure before examining outcomes.

Only 70.56% of treated participants remain at 0.10. Results therefore concern the matched treated subset, not automatically all participants. SMD < 0.10 is a diagnostic convention, not proof of exchangeability. Site-clustered intervals remain approximate and do not fully propagate propensity estimation and matching uncertainty.

The Stata teaching path uses `teffects`, a raw-scale caliper, replacement/variance conventions different from the R workflow, and an IPWRA model. Its numerical equivalence is not claimed or tested.
