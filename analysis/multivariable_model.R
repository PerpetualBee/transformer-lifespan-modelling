# =============================================================================
# multivariable_model.R
# Multi-Variable Weibull Failure Prediction Model by Cause
# Transformer Asset Lifespan Modelling Project
#
# The one-variable model (general_model.R) uses age as the sole predictor.
# This script extends the analysis by adding cause of failure as a covariate.
#
# The motivation: two transformers of the same age can have very different
# failure outlooks depending on their operating conditions. A transformer
# repeatedly struck by lightning will degrade far faster than one that ages
# naturally. Using age alone does not reveal this fully and this fuller picture is what the covariate model seeks to surface
#
# Cause of failure is treated as a categorical covariate in a Weibull
# regression. The model estimates a separate scale parameter (η) for each
# cause while sharing a single shape parameter (β) across the group.
# This is appropriate because the failure mode (how fast failure accumulates)
# varies by cause, but the general pattern of increasing risk with age (β)
# is a property of the transformer population as a whole.
#
# Representative group: 100kVA transformers
# NOTE: The same script structure was applied to the 50kVA and 200kVA groups.
#       To replicate, replace the data frame and update plot titles.
#
# Data source: Regional utility company (data sharing agreement — raw data
#              not included in this repository)
# =============================================================================


# --- 1. Setup -----------------------------------------------------------------

install.packages("survival")
install.packages("survminer")
install.packages("ggplot2")
install.packages("reshape2")

library(survival)   # Weibull MLE fitting via survreg()
library(survminer)  # Survival analysis utilities
library(ggplot2)    # Visualisation
library(reshape2)   # Data reshaping for multi-line plots


# --- 2. Data ------------------------------------------------------------------
# Same age-interval structure as the general model, with the addition of
# cause_of_failure as a categorical column. Each row's cause represents the
# dominant failure mode recorded in that age band.

HundredkVA_data <- data.frame(
  start_age        = c(0,  4,  7,   10,   13,          16,          19,                22,                  25,              28,              31),
  end_age          = c(3,  6,  9,   12,   15,          18,          21,                24,                  27,              30,              50),
  failures         = c(10, 16, 19,  10,   21,           5,           0,                 3,                   1,               1,              30),
  censored         = c(0,  0,  0,    0,    0,            0,           0,                 0,                   0,               0,               1),
  cause_of_failure = c(
    "Lightning",
    "LV Short Circuit",
    "Moisture Ingress",
    "LV Short Circuit",
    "Lightning",
    "Overloading",
    "LV Short Circuit",
    "Moisture Ingress",
    "LV Short Circuit",
    "LV Short Circuit",
    "Moisture Ingress"
  )
)


# --- 3. Format for Survival Analysis ------------------------------------------

time_to_event <- data.frame(
  age   = HundredkVA_data$end_age,
  event = ifelse(HundredkVA_data$censored == 1, 0, 1)
)

surv_obj <- Surv(time = time_to_event$age, event = time_to_event$event)


# --- 4. Fit Weibull Model with Cause of Failure as Covariate ------------------
# survreg() treats one cause as the reference level (Lightning, alphabetically
# first). The coefficients for the remaining causes are offsets from that
# reference. Exponentiating each intercept + offset gives the scale parameter
# (η) for that cause — the characteristic life specific to that failure mode.

fit_weibull_cov <- survreg(
  surv_obj ~ factor(HundredkVA_data$cause_of_failure),
  dist = "weibull"
)
summary(fit_weibull_cov)


# --- 5. Extract Parameters ----------------------------------------------------

shape <- 1 / fit_weibull_cov$scale           # β — shared across all causes
intercept <- fit_weibull_cov$coef[1]

# η for each cause: exponentiate (intercept + cause coefficient)
# Lightning is the reference level, so its η = exp(intercept) with no offset
scale_Lightning        <- exp(intercept)
scale_LV_Short_Circuit <- exp(intercept + fit_weibull_cov$coef[2])
scale_Overloading      <- exp(intercept + fit_weibull_cov$coef[3])
scale_Moisture_Ingress <- exp(intercept + fit_weibull_cov$coef[4])

cat("Shape parameter (β):", shape, "\n\n")
cat("Scale parameters (η) by cause of failure:\n")
cat("  Lightning:        ", scale_Lightning,        "years\n")
cat("  LV Short Circuit: ", scale_LV_Short_Circuit, "years\n")
cat("  Overloading:      ", scale_Overloading,      "years\n")
cat("  Moisture Ingress: ", scale_Moisture_Ingress, "years\n")

# What the scale parameters tell us:
#   A lower η means the characteristic life under that failure mode is shorter —
#   i.e. transformers exposed to that condition deteriorate faster.
#   Lightning (η ≈ 10.6 for 100kVA) versus Moisture Ingress (η ≈ 39.6) shows
#   that lightning exposure reduces expected operational life by roughly 73%.


# --- 6. Predict Failure Probabilities by Cause --------------------------------

ages <- seq(min(HundredkVA_data$start_age),
            max(HundredkVA_data$end_age), by = 1)

# Survival probability: P(transformer survives beyond age t)
# Failure probability:  1 - survival probability
survival_Lightning        <- pweibull(ages, shape = shape, scale = scale_Lightning,        lower.tail = FALSE)
survival_LV_Short_Circuit <- pweibull(ages, shape = shape, scale = scale_LV_Short_Circuit, lower.tail = FALSE)
survival_Overloading      <- pweibull(ages, shape = shape, scale = scale_Overloading,      lower.tail = FALSE)
survival_Moisture_Ingress <- pweibull(ages, shape = shape, scale = scale_Moisture_Ingress, lower.tail = FALSE)

failure_Lightning        <- 1 - survival_Lightning
failure_LV_Short_Circuit <- 1 - survival_LV_Short_Circuit
failure_Overloading      <- 1 - survival_Overloading
failure_Moisture_Ingress <- 1 - survival_Moisture_Ingress


# --- 7. Chi-Square Goodness of Fit Validation ---------------------------------
# Expected failures are calculated per cause, then summed for the overall
# expected count at each age interval. This is compared against the observed
# total using the chi-square test.

observed_failures <- HundredkVA_data$failures

expected_Lightning        <- diff(c(0, failure_Lightning))        * sum(HundredkVA_data$failures[HundredkVA_data$cause_of_failure == "Lightning"])
expected_LV_Short_Circuit <- diff(c(0, failure_LV_Short_Circuit)) * sum(HundredkVA_data$failures[HundredkVA_data$cause_of_failure == "LV Short Circuit"])
expected_Overloading      <- diff(c(0, failure_Overloading))      * sum(HundredkVA_data$failures[HundredkVA_data$cause_of_failure == "Overloading"])
expected_Moisture_Ingress <- diff(c(0, failure_Moisture_Ingress)) * sum(HundredkVA_data$failures[HundredkVA_data$cause_of_failure == "Moisture Ingress"])

expected_failures <- expected_Lightning + expected_LV_Short_Circuit +
                     expected_Overloading + expected_Moisture_Ingress

obs_exp_df <- data.frame(
  start_age = HundredkVA_data$start_age,
  end_age   = HundredkVA_data$end_age,
  observed  = observed_failures,
  expected  = tapply(
    expected_failures,
    cut(ages, breaks = c(HundredkVA_data$start_age,
                         max(HundredkVA_data$end_age))),
    sum
  )
)

chisq_test <- chisq.test(obs_exp_df$observed,
                         p         = obs_exp_df$expected,
                         rescale.p = TRUE)
print(chisq_test)


# --- 8. Plot ------------------------------------------------------------------
# Reshaping to long format allows all four cause curves to be plotted on the
# same axes with a colour legend to make the spread between causes immediately
# visible to a non-technical reader.

plot_df <- data.frame(
  age                    = ages,
  Lightning              = failure_Lightning,
  LV_Short_Circuit       = failure_LV_Short_Circuit,
  Overloading            = failure_Overloading,
  Moisture_Ingress       = failure_Moisture_Ingress
)

plot_df_melt <- reshape2::melt(plot_df, id.vars = "age")

ggplot(plot_df_melt, aes(x = age, y = value, color = variable)) +
  geom_line() +
  labs(
    title  = "Failure Probability of 100kVA Transformers by Cause",
    x      = "Age (in years)",
    y      = "Failure Probability",
    color  = "Cause of Failure"
  ) +
  scale_color_discrete(
    labels = c("Lightning", "LV Short Circuit", "Overloading", "Moisture Ingress")
  )
