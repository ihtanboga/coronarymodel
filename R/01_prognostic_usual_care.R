# coronarymodel - 3.1 Baseline usual-care prognostic design
#
# Estimand: risk stratification under usual care (NOT a treatment effect).
# Time zero = angiogram. Outcome = SCD composite. Non-SCD death is a
# competing event. Post-baseline treatment is deliberately NOT in this model.
#
# Expected columns (placeholders - swap in the real column names):
#   id
#   fu_days            days from angiography to event/censor
#   status             0 = censored, 1 = SCD/equivalent, 2 = other death
#   rca, lcx, lad, lm  0/1 territory flags
#   cad_extent         0/1/2/3 vessel disease (or richer pattern coding)
#   age, sex, diabetes, egfr, prev_mi, lvef, coronary_syndrome, year

library(data.table)
library(survival)
library(splines)

dat <- as.data.table(dat)

max_follow <- 365.25 * 10
dat[, time10   := pmin(fu_days, max_follow)]
dat[, status10 := fifelse(fu_days <= max_follow, status, 0L)]

# --- Usual-care Cox model (no post-baseline treatment) ---------------------
fit_prog <- coxph(
  Surv(time10, status10 == 1) ~
    rca + lcx + lad + lm +
    factor(cad_extent) +
    ns(age, 3) +
    sex +
    diabetes +
    ns(egfr, 3) +
    prev_mi +
    ns(lvef, 3) +
    strata(coronary_syndrome) +
    factor(year),
  data = dat, x = TRUE, y = TRUE, robust = TRUE
)

summary(fit_prog)
cox.zph(fit_prog)

# --- Standardized absolute risk (the bedside currency) ---------------------
standardized_risk <- function(fit, data, vessel, value, t = 365.25 * 5) {
  nd <- copy(data)
  nd[[vessel]] <- value
  sf <- survfit(fit, newdata = nd)
  s  <- summary(sf, times = t)$surv
  mean(1 - s, na.rm = TRUE)
}

risk_rca1_5y <- standardized_risk(fit_prog, dat, "rca", 1)
risk_rca0_5y <- standardized_risk(fit_prog, dat, "rca", 0)

risk_rca1_5y
risk_rca0_5y
risk_rca1_5y - risk_rca0_5y   # standardized 5-year absolute risk difference

# --- Cause-specific competing-risk model -----------------------------------
# install.packages("riskRegression")
library(riskRegression)

fit_csc <- CSC(
  Hist(time10, status10) ~
    rca + lcx + lad + lm +
    factor(cad_extent) +
    ns(age, 3) + sex + diabetes + ns(egfr, 3) +
    prev_mi + ns(lvef, 3) + coronary_syndrome + factor(year),
  data = dat
)

nd_rca1 <- copy(dat); nd_rca1$rca <- 1
nd_rca0 <- copy(dat); nd_rca0$rca <- 0

risk_rca1 <- predictRisk(fit_csc, newdata = nd_rca1, times = 365.25 * 5, cause = 1)
risk_rca0 <- predictRisk(fit_csc, newdata = nd_rca0, times = 365.25 * 5, cause = 1)

mean(risk_rca1)
mean(risk_rca0)
mean(risk_rca1 - risk_rca0)

# --- Split ACS from CCS -----------------------------------------------------
fit_by_syndrome <- function(syndrome_value) {
  d <- dat[coronary_syndrome == syndrome_value]
  coxph(
    Surv(time10, status10 == 1) ~
      rca + lcx + lad + lm +
      factor(cad_extent) +
      ns(age, 3) + sex + diabetes + ns(egfr, 3) +
      prev_mi + ns(lvef, 3) + factor(year),
    data = d, robust = TRUE
  )
}

fit_ccs <- fit_by_syndrome("CCS")
fit_acs <- fit_by_syndrome("ACS")

# --- ACS model with the acute variables it needs to be honest --------------
fit_acs_better <- coxph(
  Surv(time10, status10 == 1) ~
    culprit_vessel +
    rca + lcx + lad + lm +
    factor(cad_extent) +
    ns(age, 3) + sex + diabetes + ns(egfr, 3) +
    prev_mi + ns(lvef_discharge, 3) +
    peak_troponin + acute_vt_vf + reperfusion_success +
    factor(year),
  data = dat[coronary_syndrome == "ACS"], robust = TRUE
)
