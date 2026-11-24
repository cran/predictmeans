## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = TRUE, warning = TRUE)

## ----message = FALSE, warning = FALSE-----------------------------------------
library(predictmeans)

## -----------------------------------------------------------------------------
Oats <- as.data.frame(nlme::Oats)
Oats$N <- factor(Oats$nitro)
Om <- subset(Oats, !(Variety == "Victory" & N == "0.6"))
with(Om, table(Variety, N))
fm <- lmer(yield ~ N * Variety + (1 | Block/Variety), data = Om)

## ----warning = FALSE----------------------------------------------------------
pm_cells <- predictmeans(fm, "N:Variety")
nrow(pm_cells$means)

## -----------------------------------------------------------------------------
pm_n <- predictmeans(fm, "N")
pm_n$means[, c("N", "over", "estimate", "std.error", "n")]

## -----------------------------------------------------------------------------
pm_n2 <- suppressWarnings(predictmeans(fm, "N", pairwise = TRUE))
pm_n2$pairs[, c("level1", "level2", "over", "estimate", "std.error", "p.value")]

## -----------------------------------------------------------------------------
anova(fm)

## ----fig.width = 5, fig.height = 3.5------------------------------------------
ggplot2::autoplot(pm_n2, interval = "comparison")

## ----warning = FALSE----------------------------------------------------------
pm_v <- predictmeans(fm, "N:Variety", atvar = "Variety", adj = "tukey")
within <- subset(pm_v$pairs, !is.na(by))
head(within[, c("level1", "level2", "by", "estimate", "p.adj")])

## ----warning = FALSE----------------------------------------------------------
cm <- specify_contrast(Kmatrix(fm, "N:Variety"),
  N_effect_GR_and_M = (`0.6:Golden Rain` + `0.6:Marvellous`)/2 - (`0:Golden Rain` + `0:Marvellous`)/2)
contrastmeans(fm, "N:Variety", ctrmatrix = cm)[[1]]

## -----------------------------------------------------------------------------
set.seed(11)
dft <- expand.grid(rep = 1:3, trt = factor(c("A", "B")), farm = factor(c("F1", "F2")),
                   date = factor(c("d1", "d2")))
dft <- dft[!(dft$date == "d2" & dft$farm == "F2" & dft$trt == "B"), ]
dft$y <- 5 + as.integer(dft$farm) * as.integer(dft$date) + 0.5 * (dft$trt == "B") +
  rnorm(nrow(dft), 0, 0.3)
fd <- lm(y ~ date/farm + date*trt, data = dft)
predictmeans(fd, "date:trt")$means[, c("date", "trt", "estimate", "n")]

## -----------------------------------------------------------------------------
set.seed(3)
spray <- data.frame(
  droplet = factor(rep(c("medium", "fine"), c(40, 10)), levels = c("medium", "fine")),
  vol = factor(c(rep(c("130", "200"), 20), rep("400", 10))),
  conc = factor(rep(c("none", "half"), 25)))
spray$y <- 5 + 0.3 * (spray$conc == "half") - 0.4 * (spray$droplet == "fine") +
  rnorm(50, sd = 0.2)
with(spray, table(droplet, vol))

## ----warning = TRUE-----------------------------------------------------------
fs <- lm(y ~ droplet + conc + vol, data = spray)
pv <- predictmeans(fs, "vol", pairwise = TRUE)
pv$means[, c("vol", "estimate", "std.error", "n", "group")]

## -----------------------------------------------------------------------------
fs2 <- lm(y ~ conc + droplet/vol, data = spray)
pn <- predictmeans(fs2, "vol")
cbind(pn$means[, c("vol", "estimate")], fallback = pv$means$estimate,
      cells = predictmeans(fs, "vol", weights = "cells")$means$estimate)

## ----warning = TRUE-----------------------------------------------------------
set.seed(5)
past <- rbind(expand.grid(date = c("d1", "d2"), pasture = c("P1", "P2", "P3", "P4"),
                          trt = c("A", "B"), rep = 1:3),
              expand.grid(date = "d3", pasture = c("P1", "P4"), trt = c("A", "B"), rep = 1:3))
past$farm <- ifelse(past$date == "d3", ifelse(past$pasture == "P1", "F3", "F4"),
                    ifelse(past$pasture == "P4", "F2", "F1"))
for (v in c("date", "pasture", "trt", "farm")) past[[v]] <- factor(past[[v]])
past$y <- as.integer(past$date) + as.integer(past$pasture)/2 + 0.5*(past$trt == "B") +
  rnorm(nrow(past), 0, 0.3)
fp <- lm(y ~ date/farm + date/(trt + pasture)^2, data = past)
pp <- predictmeans(fp, "date:trt")

## -----------------------------------------------------------------------------
fp0 <- lm(y ~ date/(trt + pasture)^2, data = past)
p0 <- suppressWarnings(predictmeans(fp0, "date:trt"))
cbind(pp$means[, c("date", "trt", "estimate")], without_farm = p0$means$estimate)

