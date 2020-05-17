library(tidyverse)
load("genetic_application/non_gaussian_sim/combined_non_G.rda")

n <- nrow(res)
colMeans(res)
apply(res, 2, sd)/sqrt(n)
