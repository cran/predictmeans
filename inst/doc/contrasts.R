## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 4.5)
has_mvt <- requireNamespace("mvtnorm", quietly = TRUE)   # for adj = "mvt"

## -----------------------------------------------------------------------------
library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$N <- factor(Oats$nitro)
fm <- lmer(yield ~ N * Variety + (1 | Block/Variety), data = Oats)

## ----eval = has_mvt-----------------------------------------------------------
p <- sapply(c("none", "tukey", "mvt", "sidak", "scheffe", "holm"), function(a)
  predictmeans(fm, "N", adj = a)$pairs$p.adj)
pm <- predictmeans(fm, "N")
data.frame(pair = paste(pm$pairs$level1, "-", pm$pairs$level2), signif(p, 3))

## -----------------------------------------------------------------------------
pt <- predictmeans(fm, "N", adj = "tukey")$pairs
pt[c("level1", "level2", "estimate", "conf.low", "conf.high", "adj.conf.low",
     "adj.conf.high", "p.adj")]

## -----------------------------------------------------------------------------
K <- Kmatrix(fm, "N")
dunnett <- contrast_family(K, "trt.vs.ctrl")
dunnett

## ----eval = has_mvt-----------------------------------------------------------
contrastmeans(fm, "N", ctrmatrix = dunnett, adj = "mvt")

## -----------------------------------------------------------------------------
K2 <- Kmatrix(fm, "N:Variety")
trends <- contrast_family(K2, "poly", by = "Variety", max_degree = 2)
contrastmeans(fm, "N:Variety", ctrmatrix = trends)

## -----------------------------------------------------------------------------
pm_joint(fm, "N:Variety", by = "Variety")
pm_joint(fm, "N:Variety", by = "N")

## -----------------------------------------------------------------------------
curv <- contrast_family(K, "poly")[c("quadratic", "cubic"), ]
pm_joint(fm, "N", ctrmatrix = curv)

## -----------------------------------------------------------------------------
contrastmeans(fm, "N", ctrmatrix = curv, joint = TRUE)$joint

## -----------------------------------------------------------------------------
fm2 <- lmer(yield ~ nitro * Variety + (1 | Block/Variety), data = Oats)
pm_slopes(fm2, "nitro", "Variety", adj = "tukey")

## -----------------------------------------------------------------------------
pmv <- predictmeans(fm, "N:Variety")
b <- coef(pmv)
L <- (names(b) == "0.6:Victory") - (names(b) == "0.6:Golden Rain")
c(estimate = sum(L * b), SE = sqrt(drop(L %*% vcov(pmv) %*% L)))

## ----eval = requireNamespace("marginaleffects", quietly = TRUE)---------------
options(marginaleffects_model_classes = "predictmeans")
i <- match(c("0.6:Victory", "0.6:Golden Rain"), names(b))
marginaleffects::hypotheses(pmv, hypothesis = sprintf("b%d - b%d = 0", i[1], i[2]))

## -----------------------------------------------------------------------------
pe <- predictmeans(fm, "N", equivalence = 15)
pe$pairs[, c("level1", "level2", "estimate", "std.error", "p.equiv")]

## -----------------------------------------------------------------------------
fl <- lm(yield ~ Block + N + Variety, data = Oats)
pmh <- predictmeans(fl, "N", vcov = "HC3")
pmh$pairs[, c("level1", "level2", "estimate", "std.error", "p.value")]
pmh$settings$vcov
contrastmeans(fl, "N", ctrmatrix = rbind(c(-1, 1, 0, 0), c(0, -1, 1, 0)), vcov = ~ Block)

