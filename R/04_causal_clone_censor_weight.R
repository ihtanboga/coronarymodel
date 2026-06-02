# coronarymodel - 1.3 Causal design (Question B): clone-censor-weight
#
# The realistic case: no clean baseline decision, only resolved-by-90-days
# status. Clone each patient into all 3 strategies, censor a clone when it
# deviates, correct the informative censoring with IPCW. Early SCD before any
# deviation is counted in every still-adherent clone - this is what dissolves
# the "had to survive to CABG" immortal-time problem.
#
# Extra expected columns:
#   t_pci, t_cabg   days from angiography to procedure (NA if never)
#   acs_type        ACS subtype (or syndrome indicator)

library(data.table)
library(survival)
library(splines)

dat <- as.data.table(dat)

max_follow <- 365.25 * 10
dat[, time     := pmin(fu_days, max_follow)]
dat[, status10 := fifelse(fu_days <= max_follow, status, 0L)]

dat[, t_pci_i        := fifelse(is.na(t_pci), Inf, t_pci)]
dat[, t_cabg_i       := fifelse(is.na(t_cabg), Inf, t_cabg)]
dat[, t_first_revasc := pmin(t_pci_i, t_cabg_i)]

# --- Clone into the three strategies, mark deviation time ------------------
strategies <- c("PCI90", "CABG90", "MED90")
cl <- dat[rep(seq_len(nrow(dat)), each = length(strategies))]
cl[, strategy := rep(strategies, times = nrow(dat))]
cl[, dev_time := Inf]

# PCI90: first revasc within 90 days must be PCI.
cl[strategy == "PCI90" & t_pci_i  > 90, dev_time := 90]
cl[strategy == "PCI90" & t_cabg_i < pmin(t_pci_i, 90),  dev_time := pmin(dev_time, t_cabg_i)]

# CABG90: first revasc within 90 days must be CABG.
cl[strategy == "CABG90" & t_cabg_i > 90, dev_time := 90]
cl[strategy == "CABG90" & t_pci_i  < pmin(t_cabg_i, 90), dev_time := pmin(dev_time, t_pci_i)]

# MED90: no revasc within 90 days.
cl[strategy == "MED90" & t_first_revasc <= 90, dev_time := t_first_revasc]

cl[, stop        := pmin(time, dev_time)]
cl[, event_scd   := as.integer(status10 == 1 & time <= dev_time)]
cl[, cens_artif  := as.integer(dev_time < time)]

# --- IPCW for the artificial (deviation) censoring -------------------------
wdat <- survSplit(
  Surv(stop, cens_artif) ~ ., data = as.data.frame(cl),
  cut = c(30, 60, 90), start = "tstart", end = "tstop", event = "cens"
)
wdat <- as.data.table(wdat)

den_model <- glm(
  cens ~ strategy + ns(tstop, 3) +
         rca + lcx + lad + lm + factor(cad_extent) +
         ns(age, 3) + sex + diabetes + ns(egfr, 3) +
         prev_mi + ns(lvef, 3) + acs_type + factor(year),
  family = binomial(), data = wdat
)
num_model <- glm(cens ~ strategy + ns(tstop, 3), family = binomial(), data = wdat)

wdat[, p_uncens_den := 1 - predict(den_model, type = "response")]
wdat[, p_uncens_num := 1 - predict(num_model, type = "response")]
wdat[order(id, strategy, tstart, tstop),
     sw_cens := cumprod(p_uncens_num) / cumprod(p_uncens_den),
     by = .(id, strategy)]

last_w <- wdat[order(id, strategy, tstop), .SD[.N], by = .(id, strategy)][, .(id, strategy, sw_cens)]
cl <- merge(cl, last_w, by = c("id", "strategy"), all.x = TRUE)
cl[is.na(sw_cens), sw_cens := 1]
q99 <- quantile(cl$sw_cens, 0.99, na.rm = TRUE)
cl[, sw_cens_trunc := pmin(sw_cens, q99)]

# --- Weighted outcome model: interaction inside the emulated trial ---------
fit_ccw <- coxph(
  Surv(stop, event_scd) ~
    strategy * rca + strategy * lcx + strategy * lad + strategy * lm,
  data = cl, weights = sw_cens_trunc, cluster = id
)
summary(fit_ccw)

# --- Report on the absolute scale (risk-difference interaction) ------------
nd <- data.frame(
  strategy = factor(c("MED90","PCI90","CABG90","MED90","PCI90","CABG90"),
                    levels = c("MED90","PCI90","CABG90")),
  rca = c(1,1,1,0,0,0), lcx = 0, lad = 0, lm = 0
)
sf    <- survfit(fit_ccw, newdata = nd)
risk5 <- 1 - summary(sf, times = 365.25 * 5)$surv
out   <- cbind(nd, risk5)

rd_pci_vs_med_rca <- out$risk5[out$strategy=="PCI90" & out$rca==1] -
                     out$risk5[out$strategy=="MED90" & out$rca==1]
rd_pci_vs_med_no  <- out$risk5[out$strategy=="PCI90" & out$rca==0] -
                     out$risk5[out$strategy=="MED90" & out$rca==0]
interaction_rd_pci <- rd_pci_vs_med_rca - rd_pci_vs_med_no
interaction_rd_pci
