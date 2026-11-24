## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(predictmeans)
Oats <- as.data.frame(nlme::Oats)
Oats$nitro <- factor(Oats$nitro)
oats_wp <- droplevels(subset(Oats, !(Block == "I" & Variety == "Golden Rain")))
fm <- lmer(yield ~ Variety * nitro + (1 | Block/Variety), data = oats_wp)
dfs <- sapply(c("containment", "between-within", "Satterthwaite", "Kenward-Roger", "KR-expected"),
              function(m) predictmeans(fm, "Variety", df_method = m)$pairs$df)
rownames(dfs) <- c("Golden Rain - Marvellous", "Golden Rain - Victory", "Marvellous - Victory")
round(dfs, 2)

