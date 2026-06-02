# coronarymodel - 3.2 Dynamic prognostic design
#
# Estimand: residual risk RIGHT NOW given the patient's current treatment
# status. Revascularization enters as a time-varying covariate via a
# counting-process (start, stop) layout. Fixes immortal-time bookkeeping;
# does NOT remove confounding by indication. Prognostically valid,
# causally mute.
#
# Extra expected columns:
#   scd_event   0/1 SCD indicator at fu_days
#   pci_day     day of PCI  (NA if none)
#   cabg_day    day of CABG (NA if none)

library(survival)
library(splines)

# tmerge builds the (tstart, tstop) structure and flips treatment status at
# the actual procedure date.
td <- tmerge(
  data1 = dat, data2 = dat, id = id,
  scd = event(fu_days, scd_event)
)

td <- tmerge(
  td, dat, id = id,
  pci_td  = tdc(pci_day),
  cabg_td = tdc(cabg_day)
)

td$tx_td <- with(td,
  ifelse(cabg_td == 1, "CABG",
  ifelse(pci_td  == 1, "PCI", "NO_REVASC"))
)
td$tx_td <- factor(td$tx_td, levels = c("NO_REVASC", "PCI", "CABG"))

fit_dynamic <- coxph(
  Surv(tstart, tstop, scd) ~
    rca + lcx + lad + lm +
    tx_td +
    factor(cad_extent) +
    ns(age, 3) + sex + diabetes + ns(egfr, 3) +
    prev_mi + ns(lvef, 3) +
    strata(coronary_syndrome),
  data = td, cluster = id
)

summary(fit_dynamic)
