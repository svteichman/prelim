library(tidyverse)
load("simulations/perm_sim/permutation_res.rda")

n <- dim(perm_sim)[2]
mean(unlist(perm_sim[1,1,]))
mean(unlist(perm_sim[1,2,]))
mean(unlist(perm_sim[1,3,]))
mean(unlist(perm_sim[1,4,]))
mean(unlist(perm_sim[1,5,]))
sd(unlist(perm_sim[1,1,]))/sqrt(n)
sd(unlist(perm_sim[1,2,]))/sqrt(n)
sd(unlist(perm_sim[1,3,]))/sqrt(n)
sd(unlist(perm_sim[1,4,]))/sqrt(n)
sd(unlist(perm_sim[1,5,]))/sqrt(n)

load("simulations/perm_sim/permutation_res_plus.rda")

n <- dim(perm_sim1)[2]
mean(unlist(perm_sim1[1,1,]))
mean(unlist(perm_sim1[1,2,]))
mean(unlist(perm_sim1[1,3,]))
mean(unlist(perm_sim1[1,4,]))
mean(unlist(perm_sim1[1,5,]))
sd(unlist(perm_sim1[1,1,]))/sqrt(n)
sd(unlist(perm_sim1[1,2,]))/sqrt(n)
sd(unlist(perm_sim1[1,3,]))/sqrt(n)
sd(unlist(perm_sim1[1,4,]))/sqrt(n)
sd(unlist(perm_sim1[1,5,]))/sqrt(n)
