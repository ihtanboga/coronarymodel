# coronarymodel - 1.3 Causal design (Question B): IPTW ITT
#
# Use ONLY when a genuine baseline intended strategy (tx0) exists in the data,
# decided at angiography - NOT reconstructed from the three-month window.
# Multinomial propensity, time-zero confounders only, stabilized truncated
# weights. The tx0 x vessel interactions formally test effect modification.
#
# Requires time10 / status10 (see 01_prognostic_usual_care.R for their
# construction from fu_days / status).
#
# Extra expected columns:
#   significant_cad   0/1
#   tx0               baseline intended strategy: MEDICAL / PCI / CABG
#   acs_type          ACS subtype (or syndrome indicator)

library(nnet)
library(survival)
library(splines)
library(data.table)

trial <- copy(dat)
trial <- trial[significant_cad == 1 & !is.na(tx0)]
trial[, tx0 := factor(tx0, levels = c("MEDICAL", "PCI", "CABG"))]

# Multinomial propensity for treatment assignment, time-zero covariates only.
ps_model <- multinom(
  tx0 ~ rca + lcx + lad + lm + factor(cad_extent) +
        ns(age, 3) + sex + diabetes + ns(egfr, 3) +
        prev_mi + ns(lvef, 3) + acs_type + factor(year),
  data = trial, trace = FALSE
)

ps_denom <- fitted(ps_model)
idx      <- cbind(seq_len(nrow(trial)), match(trial$tx0, colnames(ps_denom)))
p_denom  <- ps_denom[idx]
p_num    <- as.numeric(prop.table(table(trial$tx0))[as.character(trial$tx0)])

trial[, sw_tx := p_num / p_denom]
q99 <- quantile(trial$sw_tx, 0.99, na.rm = TRUE)
trial[, sw_tx_trunc := pmin(sw_tx, q99)]   # truncation matters in practice

fit_itt <- coxph(
  Surv(time10, status10 == 1) ~
    tx0 * rca + tx0 * lcx + tx0 * lad + tx0 * lm,
  data = trial, weights = sw_tx_trunc, cluster = id
)
summary(fit_itt)
