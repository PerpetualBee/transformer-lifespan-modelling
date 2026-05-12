# =============================================================================
# general_model.R
# One-Variable Weibull Failure Prediction Model
# Transformer Asset Lifespan Modelling Project
#
# This script fits a 2-parameter Weibull distribution to transformer failure
# data using age as the sole predictor. It estimates shape (β) and scale (η)
# parameters, predicts failure probabilities over the transformer lifecycle,
# and validates the model with a Chi-Square Goodness of Fit Test.
#
# Representative group: 100kVA transformers
# NOTE: The same script structure was applied to the 50kVA and 200kVA groups.
#       To replicate for another group, replace the data frame below and
#       update the plot title accordingly.
#
# Data source: Regional utility company (data sharing agreement; raw data
#              not included in this repository)
# =============================================================================


# --- 1. Setup -----------------------------------------------------------------

install.packages("survival")
install.packages("survminer")
install.packages("ggplot2")

library(survival)   # Weibull MLE fitting via survreg()
library(survminer)  # Survival analysis utilities
library(ggplot2)    # Visualisation


# --- 2. Data ------------------------------------------------------------------
# Failure records are organised into age intervals with recorded failure counts.
# The final interval is right-censored (censored = 1), meaning some transformers
# in that range had not yet failed by the end of the observation window.
# Treating these as censored rather than failures avoids underestimating
# the characteristic life.

HundredkVA_Dataset <- data.frame(
  Start_Age         = c(0,  4,  7,  10, 13, 16, 19, 22, 25, 28, 31),
  End_Age           = c(3,  6,  9,  12, 15, 18, 21, 24, 27, 30, 50),
  Recorded_Failures = c(10, 16, 19, 10, 21,  5,  0,  3,  1,  2, 30),
  Censored          = c(0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  1)
)


# --- 3. Format for Survival Analysis ------------------------------------------
# survreg() requires a single time-to-event format.
# Each row becomes one observation: the end age of the interval and whether
# a failure (event = 1) or censoring (event = 0) occurred.

time_to_event <- data.frame(
  age   = HundredkVA_Dataset$End_Age,
  event = ifelse(HundredkVA_Dataset$Censored == 1, 0, 1)
)

surv_obj <- Surv(time = time_to_event$age, event = time_to_event$event)


# --- 4. Fit Weibull Distribution ----------------------------------------------
# The Weibull distribution was chosen over normal and lognormal alternatives
# because it is the industry-standard method for modelling asset failure data.
# It is flexible enough to capture all three stages of the classic bathtub
# curve depending on the value of β, making it well-suited to transformer
# lifetime data where the failure mode is not assumed in advance.
#
# survreg() uses Maximum Likelihood Estimation (MLE) to find the parameters
# that make the observed failure data most probable given the Weibull model.
# Note: survreg() returns parameters in a different parameterisation, so
# shape and scale are extracted with the transformations below.

fit_weibull <- survreg(surv_obj ~ 1, dist = "weibull")
summary(fit_weibull)

shape <- 1 / fit_weibull$scale     # β: shape parameter
scale <- exp(fit_weibull$coef)     # η: scale parameter (characteristic life)

cat("Shape parameter (β):", shape, "\n")
cat("Scale parameter (η):", scale, "years\n")

# Interpretation guide:
#   β > 1  →  failure rate is increasing (wear-out behaviour)
#   β = 1  →  constant failure rate (random failures only)
#   β < 1  →  decreasing failure rate (early-life / infant mortality)
#
#   η      →  the age at which ~63.2% of the population has failed.
#             In a maintenance context, this is the age at which a sharp
#             uptick in failures begins, basically the trigger point for
#             proactive replacement programmes.


# --- 5. Predict Failure Probabilities -----------------------------------------

ages           <- seq(min(HundredkVA_Dataset$Start_Age),
                      max(HundredkVA_Dataset$End_Age), by = 1)
survival_probs <- pweibull(ages, shape = shape, scale = scale, lower.tail = FALSE)
failure_probs  <- 1 - survival_probs


# --- 6. Chi-Square Goodness of Fit Validation ---------------------------------
# Compares observed failure counts with the counts the model predicts at each
# age interval. A result above the critical value (18.307 at α = 0.05, df = 10)
# means the model's predicted distribution differs significantly from the
# observed data. This does not mean Weibull is wrong for this data as it
# remained the best fit compared to normal and lognormal alternatives (this was tested as well).
# The deviation is attributed to the short observation window (30 months)
# and the absence of diagnostic covariates in the dataset.

observed_failures <- HundredkVA_Dataset$Recorded_Failures
expected_failures <- diff(c(0, failure_probs)) * sum(HundredkVA_Dataset$Recorded_Failures)

obs_exp_df <- data.frame(
  start_age = HundredkVA_Dataset$Start_Age,
  end_age   = HundredkVA_Dataset$End_Age,
  observed  = observed_failures,
  expected  = tapply(
    expected_failures,
    cut(ages, breaks = c(HundredkVA_Dataset$Start_Age,
                         max(HundredkVA_Dataset$End_Age))),
    sum
  )
)

chisq_test <- chisq.test(obs_exp_df$observed,
                         p          = obs_exp_df$expected,
                         rescale.p  = TRUE)
print(chisq_test)


# --- 7. Plot ------------------------------------------------------------------

plot_df <- data.frame(age = ages, failure_prob = failure_probs)

ggplot(plot_df, aes(x = age, y = failure_prob)) +
  geom_line() +
  labs(
    title = "Failure Probability of 100kVA Power Transformers",
    x     = "Age (in years)",
    y     = "Failure Probability"
  )
