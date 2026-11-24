## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 6)

## -----------------------------------------------------------------------------
library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
fm <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
pm_check(fm)

## ----fig.alt = "Four diagnostic plots of the split-plot model"----------------
residplot(fm)

## -----------------------------------------------------------------------------
gm <- glm(breaks ~ wool * tension, family = poisson, data = warpbreaks)
pm_check(gm)

## ----fig.alt = "Diagnostic plots of the Poisson model"------------------------
residplot(gm)

## -----------------------------------------------------------------------------
nb <- MASS::glm.nb(breaks ~ wool * tension, data = warpbreaks)
pm_check(nb)

## -----------------------------------------------------------------------------
pm_nb <- predictmeans(nb, "wool:tension", atvar = "tension", pairwise = TRUE)
pm_nb$pairs[!is.na(pm_nb$pairs$by), c("level1", "level2", "ratio", "ratio.conf.low",
                                      "ratio.conf.high", "p.value")]

## ----eval = FALSE-------------------------------------------------------------
# residplot(fm, interactive = TRUE)

