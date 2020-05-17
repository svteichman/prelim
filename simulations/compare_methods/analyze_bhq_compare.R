library(tidyverse)
load("simulations/compare_methods/combined_bhq.rda")

n <- nrow(res)
colMeans(res)
apply(res, 2, sd)/sqrt(n)
