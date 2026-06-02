# coronarymodel - 1.3 Landmark sensitivity + completeness of revascularization
#
# A 90-day landmark is a fair SENSITIVITY check only: it classifies treatment
# among those alive and event-free at day 90, so it deletes the early deaths
# that matter most. Report those early events separately; this is never the
# primary causal answer.
#
# Then a reframing closer to the real clinical question: completeness of
# revascularization (were ALL significant vessels actually fixed), which is the
# more honest exposure if the RCA signal is really residual disease.
#
# Requires time / status10 (see 04_causal_clone_censor_weight.R for their
# construction from fu_days / status).
#
# Extra expected columns:
#   rca_treated90, lcx_treated90, lad_treated90, lm_treated90  0/1 by day 90

library(data.table)
library(survival)
library(splines)

dat <- as.data.table(dat)
dat[, t_pci_i  := fifelse(is.na(t_pci),  Inf, t_pci)]
dat[, t_cabg_i := fifelse(is.na(t_cabg), Inf, t_cabg)]

# --- 90-day landmark --------------------------------------------------------
land <- copy(dat)[time > 90]
land[, t_land := time - 90]
land[, tx90 := fifelse(t_cabg_i <= 90, "CABG",
              fifelse(t_pci_i  <= 90, "PCI", "MEDICAL"))]
land[, tx90 := factor(tx90, levels = c("MEDICAL", "PCI", "CABG"))]

fit_landmark <- coxph(
  Surv(t_land, status10 == 1) ~
    tx90 * rca + tx90 * lcx + tx90 * lad + tx90 * lm +
    factor(cad_extent) + ns(age, 3) + sex + diabetes + ns(egfr, 3) +
    prev_mi + ns(lvef, 3) + strata(coronary_syndrome),
  data = land, cluster = id
)
summary(fit_landmark)

# --- Completeness of revascularization by day 90 ---------------------------
# Complete = every significant vessel was treated.
dat[, complete90 := as.integer(
  (rca == 0 | rca_treated90 == 1) &
  (lcx == 0 | lcx_treated90 == 1) &
  (lad == 0 | lad_treated90 == 1) &
  (lm  == 0 | lm_treated90  == 1)
)]
