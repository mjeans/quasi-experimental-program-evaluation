version 18
clear all
set more off

* Companion implementation for the R-generated synthetic dataset.
* Run R/01_generate_synthetic_data.R before this do-file.

capture confirm file "data/synthetic_program_data.csv"
if _rc {
    display as error "Generate data/synthetic_program_data.csv first."
    exit 601
}

import delimited using "data/synthetic_program_data.csv", clear varnames(1)
encode site_id, gen(site_n)
encode region, gen(region_n)

isid participant_id
assert inlist(program, 0, 1)
misstable summarize outcome baseline_score age female urban prior_service

* Naive comparison, shown as a benchmark rather than a causal estimate.
regress outcome i.program, vce(cluster site_n)

* Propensity-score matching for the average effect on participants.
teffects psmatch (outcome) ///
    (program c.baseline_score c.age i.female i.urban ///
        i.prior_service i.region_n), ///
    atet nneighbor(1) caliper(.20) vce(robust)
estimates store psm_att
tebalance summarize
tebalance density baseline_score

* Doubly robust ATT estimate as a specification check.
teffects ipwra ///
    (outcome c.baseline_score c.age i.female i.urban ///
        i.prior_service i.region_n) ///
    (program c.baseline_score c.age i.female i.urban ///
        i.prior_service i.region_n), ///
    atet vce(robust)
estimates store ipwra_att

estimates table psm_att ipwra_att, b(%7.3f) se(%7.3f) stats(N)
