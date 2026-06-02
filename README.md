# coronarymodel

This repo is the methods companion to a re-analysis argument: when a coronary
cohort starts the clock at angiography and then reads a three-month treatment
status as if it modified vessel-specific sudden-death risk, the model is
answering a different question than the one being claimed. The fix is not a
better single Cox model. The fix is to decide which question you are answering
and then build the design for *that* question.

Below I lay out the prognostic design, the dynamic-prediction design, and the
causal design. Each one is implemented as a runnable R skeleton under
[`R/`](./R); the scripts are referenced from the relevant section here.
Variable names are placeholders; swap in the real column names from your data.

A note on what is shared across all of them: time zero is the angiogram, the
outcome is the SCD composite (ideally with the components also reported
separately), and non-SCD death is a competing event, not a censoring event.

---

## Part III: What a Correct Re-Analysis Would Look Like

### 3.1 The prognostic design (Question A), done properly

This is the section the original draft treated too thinly, so here it is in
full.

The goal is risk stratification under usual care, not a treatment effect. That
single decision dictates the whole structure, and in particular it dictates
that the three-month treatment variable does **not** go into the main vessel
model as a confounder. Treatment is downstream of the anatomy. Putting it in
answers a different question (the conditioned residual one) and pulls in the
collider problem from Part II.3.

The protocol:

- **Population:** first angiogram, significant CAD present.
- **Time zero:** date of first angiogram.
- **Exposure:** the four territory flags, plus an explicit extent variable
  (1/2/3-vessel, LM±others), because the supplement told us extent is where the
  real signal lives.
- **Outcome:** fatal SCD, aborted SCD, and ICD therapy, reported separately and
  as the composite.
- **Covariates:** baseline-only. Age, sex, diabetes, eGFR, prior MI, syndrome
  type, baseline LVEF, treatment year. No post-baseline treatment.
- **Reporting:** hazard ratios are not enough. Report standardized 5- and
  10-year cumulative incidence, and do it under a competing-risk framework.
- **ACS and CCS analyzed separately.**

The implementation is in
[`R/01_prognostic_usual_care.R`](R/01_prognostic_usual_care.R). It fits a
usual-care Cox model with baseline-only covariates and strata for syndrome
type, deliberately excluding the post-baseline treatment status. That model
answers: within the same baseline profile and disease extent, is an
RCA/LCx/LAD/LM flag associated with SCD?

But a hazard ratio is still the wrong currency for a bedside decision. The same
script converts the fit to standardized absolute risk, so that the sentence
becomes "RCA stenosis raises standardized 5-year SCD risk by X percentage
points," which a clinician can actually use, instead of "HR 1.53," which they
cannot.

When non-SCD death is common (and in a cohort this old and this sick, it is),
Cox alone understates the competing-risk reality, so the script also fits a
cause-specific competing-risk model and reports standardized 5-year risk by
vessel flag on the absolute scale.

Finally it splits ACS from CCS rather than pooling them. For ACS the model is
only honest if the acute variables are in it — culprit status, discharge LVEF,
troponin, acute arrhythmia, reperfusion success. Without those, an ACS vessel
coefficient is close to uninterpretable.

---

### 3.2 The dynamic prognostic design

The user specifically asked for this, and it is worth separating from both the
static prognostic model and the causal model, because it answers yet another
distinct question: *as time passes and the patient's treatment status updates,
what is their residual risk right now?*

This is the legitimate way to let revascularization into a prognostic model.
You do not pretend to know at angiography whether someone will be
revascularized. You let the covariate change value when the procedure actually
happens, using a counting-process (start, stop) data layout. The implementation
is in [`R/02_dynamic_prediction.R`](R/02_dynamic_prediction.R), which uses
`tmerge` to build the `(tstart, tstop)` structure and flip treatment status at
the actual procedure date.

The output is useful for live prediction: "given that this patient has now had
PCI, here is the residual SCD hazard." But, and this is the boundary the
original draft did not draw clearly enough, the time-varying treatment
coefficient is still **not** a causal treatment effect. Treatment may have been
chosen *because* the patient deteriorated, *because* the anatomy turned out
complex, *because* the physician saw trouble coming. The time-varying covariate
fixes the immortal-time bookkeeping (because status only flips when the
procedure truly happens, nobody is credited with treatment time they did not
live through), but it does nothing about confounding by indication.
Prognostically valid, causally mute.

This is also, incidentally, what the paper's own Supplementary Table 6 gestures
at, except that table still keeps the three-month decision covariate in the
model, which reintroduces exactly the conditioning we are trying to avoid.

---

### 3.3 The causal design (Question B): target-trial emulation

If you want to claim that treatment strategy modifies vessel-specific risk,
this is the bar you have to clear. Not a stratified Cox. A target trial.

The protocol, written as a trial first and then emulated:

| Component | Specification |
| --- | --- |
| Eligibility | Significant CAD at first angiogram; candidate for at least one strategy |
| Time zero | Angiography / heart-team decision date |
| Strategies | (1) PCI within 90 days; (2) CABG within 90 days; (3) no revascularization within 90 days / OMT |
| Assignment | Randomization in a real trial; in emulation, baseline decision or clone-censor-weight |
| Follow-up | From time zero, not from the date the procedure happened |
| Outcome | SCD, fatal SCD, aborted SCD, ICD therapy; competing death handled separately |
| Contrast | Standardized 5-year absolute risk difference and risk ratio; formal strategy × vessel interaction |
| Confounders | Time-zero variables only |
| Analysis | IPTW / g-formula / MSM; clone-censor-weight if a grace period exists |

**If a genuine baseline intended strategy exists** in the data (decided at
angiography, not reconstructed from the three-month window), an ITT-style
weighted analysis is reasonable. This is implemented in
[`R/03_causal_iptw_itt.R`](R/03_causal_iptw_itt.R): a multinomial propensity
model on time-zero covariates only, stabilized and truncated weights, and an
outcome model whose `strategy × vessel` interaction terms formally test "does
treatment effect depend on vessel" — the question Figure 1 only pretends to
answer. Even here the causal reading requires that all the important baseline
confounders were actually measured, which, given the missing drug data, is not
obviously true.

**The more realistic situation** in this dataset is that there is no clean
baseline decision, only the resolved-by-90-days status. In that case the
appropriate primary analysis is clone-censor-weight, which is the proper way to
handle a grace period without manufacturing immortal time. Each patient is
cloned into all three strategies, each clone is censored when it deviates from
its assigned strategy, and the informative censoring is corrected with
inverse-probability-of-censoring weights. This is implemented in
[`R/04_causal_clone_censor_weight.R`](R/04_causal_clone_censor_weight.R).

The point of cloning is that an early SCD, occurring before any deviation, gets
counted in every clone that has not yet deviated. That is precisely what
dissolves the "you had to survive to CABG to be in the CABG group" problem. The
script then fits the weighted outcome model — testing the interaction inside an
emulated trial, not across selected strata — and reports it on the absolute
scale, because an interaction that looks absent on the hazard-ratio scale can be
large and clinically real on the risk-difference scale.

A 90-day landmark is a fair sensitivity check (it cuts the post-baseline
misclassification by only classifying treatment among those alive and
event-free at day 90), but it cannot be the primary causal answer, because it
deletes the early deaths that matter most. Report those early events
separately. The landmark sensitivity model lives in
[`R/05_landmark_and_completeness.R`](R/05_landmark_and_completeness.R).

One more reframing that I think is closer to the real clinical question than
the PCI/CABG/medical label: completeness of revascularization. If the RCA
signal is really a residual-disease signal, then the variable that matters is
not "which procedure" but "were all the significant vessels actually fixed."
That, too, is post-baseline and needs a landmark or clone-censor-weight
treatment, but it is the more honest exposure. The completeness flag is
constructed at the end of the same landmark script.

---

### 3.4 The short version of the design argument

The "should revascularization be in the model" question has no single answer.
It depends entirely on the estimand:

| Goal | Revascularization in the main model? | Correct handling |
| --- | --- | --- |
| Baseline usual-care prognosis | No | Report treatment as a downstream process |
| Baseline conditional prognosis | Only if the baseline strategy is truly known at baseline | "Risk by planned treatment" |
| 90-day post-treatment prognosis | Yes | Landmark model |
| Dynamic prediction | Yes | Time-varying covariate |
| Causal PCI/CABG/medical comparison | Yes | Target-trial emulation / IPTW / clone-censor-weight |
| "Is the RCA biologically more dangerous?" | Generally no | Treatment is mediator and collider |

Rankinen's choice, putting a three-months-resolved treatment status into a Cox
model whose clock started at angiography, and then reading the strata as causal
modification, sits in none of these cells cleanly. It mixes the bookkeeping of
one estimand with the interpretation of another.

---

## Scripts

| File | Estimand |
| --- | --- |
| [`R/01_prognostic_usual_care.R`](R/01_prognostic_usual_care.R) | Baseline usual-care prognosis (3.1) |
| [`R/02_dynamic_prediction.R`](R/02_dynamic_prediction.R) | Dynamic / residual risk, time-varying treatment (3.2) |
| [`R/03_causal_iptw_itt.R`](R/03_causal_iptw_itt.R) | Causal ITT with baseline strategy, IPTW (3.3) |
| [`R/04_causal_clone_censor_weight.R`](R/04_causal_clone_censor_weight.R) | Causal with grace period, clone-censor-weight (3.3) |
| [`R/05_landmark_and_completeness.R`](R/05_landmark_and_completeness.R) | 90-day landmark sensitivity + completeness of revascularization (3.3) |

All variable names are placeholders. Swap in the real column names from your
data before running.

### Dependencies

```
install.packages(c("data.table", "survival", "splines",
                   "riskRegression", "nnet"))
```
