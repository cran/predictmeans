## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(predictmeans)
args(R2_model)

## -----------------------------------------------------------------------------
data(Body_length)
str(Body_length)
mod_BodyL <- lmer(BodyL ~ Treatment + Habitat + Sex + (1|Population/Container),
                  data = Body_length)
R2_model(mod_BodyL)

## -----------------------------------------------------------------------------
data(Poisson_ind)
mod_poisson <- glmer(Egg ~ Treatment + Habitat + (1|Population/Container/individual),
                     family = "poisson", data = Poisson_ind)
R2_model(mod_poisson, over_disp = TRUE)

## -----------------------------------------------------------------------------
data(Clinical)
Clinical$Unfavorable <- Clinical$Total - Clinical$Favorable
fit_grp <- function(link)
  glmer(cbind(Favorable, Unfavorable) ~ Treatment + (1|Clinic/Treatment),
        family = binomial(link = link), data = Clinical)
R2_model(fit_grp("logit"))
R2_model(fit_grp("probit"))
R2_model(fit_grp("cloglog"))

## -----------------------------------------------------------------------------
idx <- rep(seq_len(nrow(Clinical)), Clinical$Total)
Clinical_N <- Clinical[idx, c("Clinic", "Treatment")]
Clinical_N$favour <- unlist(lapply(seq_len(nrow(Clinical)), function(i)
  rep(c(1, 0), c(Clinical$Favorable[i], Clinical$Unfavorable[i]))))
str(Clinical_N)
fit_pat <- function(link)
  glmer(favour ~ Treatment + (1|Clinic/Treatment),
        family = binomial(link = link), data = Clinical_N)
R2_model(fit_pat("logit"))
R2_model(fit_pat("probit"))
R2_model(fit_pat("cloglog"))

