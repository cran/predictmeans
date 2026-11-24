## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
try(predictmeans(fm, "nitro", level = 0.1, Df = 30))

## ----eval = FALSE-------------------------------------------------------------
# fm0 <- lmer(yield ~ nitro + Variety + (1 | Block/Variety), data = Oats)
# fmc <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = nlme::Oats)
# fs <- semireg(rate ~ s(pressure, by = QB, bs = "os", k = 4) + (pressure | Subject),
#               data = nlme::Dialyzer)

## -----------------------------------------------------------------------------
Om <- Oats[-c(1, 2, 30), ]   # unbalanced
fm_u <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Om)
new <- predictmeans(fm_u, "Variety")
old <- predictmeans(fm_u, "Variety", vcov_adjust = FALSE)
cbind(SE_2.0 = new$means$std.error, SE_1.x = old$means$std.error)

## -----------------------------------------------------------------------------
pm <- predictmeans(fm, "nitro", adj = "BH")
names(pm)
pm$`Predicted Means`
pm$means

