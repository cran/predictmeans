## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 4)

## -----------------------------------------------------------------------------
library(predictmeans)

## -----------------------------------------------------------------------------
data(Drug)
fm_drug <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
pm_drug <- predictmeans(fm_drug, "Drug", pairwise = TRUE)
pm_drug$means[, c("Drug", "estimate", "response", "response.std.error",
                  "response.conf.low", "response.conf.high")]

## -----------------------------------------------------------------------------
pm_drug$pairs[, c("level1", "level2", "estimate", "p.value",
                  "odds.ratio", "odds.ratio.conf.low", "odds.ratio.conf.high")]

## -----------------------------------------------------------------------------
ggplot2::autoplot(pm_drug, type = "ci")

## -----------------------------------------------------------------------------
data(Poisson_ind)
fm_egg <- glmer(Egg ~ Treatment + Habitat + (1 | Population/Container/individual),
                family = poisson, data = Poisson_ind)
pm_egg <- predictmeans(fm_egg, "Treatment", pairwise = TRUE)
pm_egg$means[, c("Treatment", "estimate", "response", "response.std.error",
                 "response.conf.low", "response.conf.high")]
pm_egg$pairs[, c("level1", "level2", "p.value", "ratio", "ratio.conf.low", "ratio.conf.high")]

## -----------------------------------------------------------------------------
pm_pop <- predictmeans(fm_egg, "Treatment", bias_adjust = TRUE, average = "response")
pm_pop$means[, c("Treatment", "response", "response.std.error")]
pm_pop$settings[c("bias_adjust", "average")]

## -----------------------------------------------------------------------------
data(Body_length)
fm_len <- lmer(log(BodyL) ~ Treatment + Habitat + Sex + (1 | Population/Container),
               data = Body_length)
pm_len <- predictmeans(fm_len, "Treatment", pairwise = TRUE, trans = exp)
pm_len$means[, c("Treatment", "estimate", "response", "response.std.error")]
pm_len$pairs[, c("level1", "level2", "p.value", "ratio", "ratio.conf.low", "ratio.conf.high")]

## -----------------------------------------------------------------------------
data(Clinical)
fm_clin <- glmer(cbind(Favorable, Total - Favorable) ~ Treatment + (1 | Clinic),
                 family = binomial, data = Clinical)
pm_clin <- predictmeans(fm_clin, "Treatment", pairwise = TRUE)
pm_clin$means[, c("Treatment", "response", "response.conf.low", "response.conf.high")]
pm_clin$pairs[, c("level1", "level2", "p.value", "odds.ratio",
                  "odds.ratio.conf.low", "odds.ratio.conf.high")]

## -----------------------------------------------------------------------------
R2_model(fm_clin)
R2_model(fm_clin, method = "nakagawa")

## ----fig.width = 6, fig.height = 4--------------------------------------------
set.seed(1)
dose <- data.frame(trt = gl(2, 40, labels = c("control", "treated")), x = runif(80, 0, 2))
dose$count <- rpois(80, exp(1 + 0.2*dose$x + 0.6*dose$x*(dose$trt == "treated")))
fm_dose <- glm(count ~ trt*x, family = poisson, data = dose)
cmp <- covariatemeans(fm_dose, "trt", covariate = "x", contr = c("treated", "control"),
                      band = "simultaneous")
cmp
ggplot2::autoplot(cmp)

