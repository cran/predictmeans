## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 4.5)

## ----echo = FALSE, out.width = "100%", fig.alt = "From a fitted model to the report: diagnostics, tests, means and contrasts, variance components, then tables and plots"----
knitr::include_graphics("../man/figures/workflow.png")

## -----------------------------------------------------------------------------
library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)

## -----------------------------------------------------------------------------
pm <- predictmeans(fm, "nitro:Variety", atvar = "Variety", adj = "tukey")
summary(pm)

## -----------------------------------------------------------------------------
head(pm$means)
head(pm$pairs)

## -----------------------------------------------------------------------------
head(as.data.frame(pm, what = "pairs"), 3)

## -----------------------------------------------------------------------------
ggplot2::autoplot(pm)
ggplot2::autoplot(pm, type = "pvalue")

## -----------------------------------------------------------------------------
ggplot2::autoplot(pm, type = "bar")
ggplot2::autoplot(pm, type = "ci", data = TRUE)

## -----------------------------------------------------------------------------
# Satterthwaite df
pm_s <- predictmeans(fm, "Variety", df_method = "Satterthwaite")
# SAS DDFM=CONTAIN (containment) df, with unadjusted standard errors
pm_c <- predictmeans(fm, "Variety", df_method = "containment")
rbind(Satterthwaite = pm_s$means$df[1], containment = pm_c$means$df[1])

## -----------------------------------------------------------------------------
predictmeans(fm, "nitro", weights = "proportional")$means

## -----------------------------------------------------------------------------
cm <- specify_contrast(Kmatrix(fm, "nitro:Variety"),
  N_effect_Victory = `0.6:Victory` - `0:Victory`,
  N_effect_Victory_vs_others = (`0.6:Victory` - `0:Victory`) -
    ((`0.6:Golden Rain` - `0:Golden Rain`) + (`0.6:Marvellous` - `0:Marvellous`)) / 2)
contrastmeans(fm, "nitro:Variety", ctrmatrix = cm)

## -----------------------------------------------------------------------------
data(Drug)
fm_glm <- glm(cbind(R, N - R) ~ Drug + log10Dose, family = binomial, data = Drug)
predictmeans(fm_glm, "Drug", pairwise = TRUE)$means

