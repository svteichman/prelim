library(prelim.knockoffs)
library(mvtnorm)
library(knockoff)

# Knockoff 

normalize <- function(X) {
  X.scaled <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  X_scaled <- X.scaled[, ]
  return(X_scaled)
}

generate_dat <- function() {
  n <- 300
  p <- 100
  mean <- rep(0,p)
  cov <- matrix(data = 0.3, nrow = p, ncol = p) + diag(x = 0.7, nrow = p)
  X <- rmvnorm(n, mean, cov)
  X <- normalize(scale(X, center = T, scale = F))
  z <- rnorm(n, 0, 1) 
  x_30_sum <- apply(X, 1, function(x) sum(x[1:30]))
  y <- 3.5*x_30_sum + z
  return(list(X,y))
}

get_fdr <- function(plus) {
  dat <- generate_dat()
  X <- dat[[1]]
  y <- dat[[2]]
  knock_sdp_selected <- knockoff_filter(X, y, "sdp", fdr = .2, plus = plus, randomize = TRUE)$selected
  knock_equi_selected <- knockoff_filter(X, y, "equi", fdr = .2, plus = plus, randomize = TRUE)$selected
  
  offset <- as.numeric(plus)
  
  knock.gen.equi <- function(x) create.fixed(x, method = "equi", randomize = TRUE)
  knock.gen.sdp <- function(x) create.fixed(x, method = "sdp", randomize = TRUE)
  their_e <- knockoff.filter(X, y, fdr = .2, knockoffs=knock.gen.equi, 
                             statistic=stat.glmnet_lambdasmax, offset = offset)$selected
  their_s <- knockoff.filter(X, y, fdr = .2, knockoffs=knock.gen.sdp, 
                             statistic=stat.glmnet_lambdasmax, offset = offset)$selected
  
  perm_selected <- permutation_filter(X, y, fdr = .2, plus = plus)$selected
  
  knock_sdp_fdr <- sum(knock_sdp_selected > 30)/(max(1,length(knock_sdp_selected)))
  knock_equi_fdr <- sum(knock_equi_selected > 30)/(max(1,length(knock_equi_selected)))
  their_s_fdr <- sum(their_s > 30)/(max(1,length(their_s)))
  their_e_fdr <- sum(their_e > 30)/(max(1,length(their_e)))
  perm_fdr <- sum(perm_selected > 30)/(max(1,length(perm_selected)))
  res <- data.frame(my_s_k = knock_sdp_fdr, my_e_k = knock_equi_fdr, their_s_k = their_s_fdr,
                    their_e_k = their_e_fdr, perm = perm_fdr)
  return(res)
}
n <- 1000
sim <- replicate(n, get_fdr(plus = FALSE))
mean(unlist(sim[1,1,]))
mean(unlist(sim[1,2,]))
mean(unlist(sim[1,3,]))
mean(unlist(sim[1,4,]))
mean(unlist(sim[1,5,]))
sd(unlist(sim[1,1,]))/sqrt(n)
sd(unlist(sim[1,2,]))/sqrt(n)
sd(unlist(sim[1,3,]))/sqrt(n)
sd(unlist(sim[1,4,]))/sqrt(n)
sd(unlist(sim[1,5,]))/sqrt(n)

perm_sim <- sim
save(perm_sim, file = "simulations/perm_sim/permutation_res.rda")


# trying with Knockoff+ 

sim1 <- replicate(n, get_fdr(plus = TRUE))
mean(unlist(sim1[1,1,]))
mean(unlist(sim1[1,2,]))
mean(unlist(sim1[1,3,]))
mean(unlist(sim1[1,4,]))
mean(unlist(sim1[1,5,]))
sd(unlist(sim1[1,1,]))/sqrt(n)
sd(unlist(sim1[1,2,]))/sqrt(n)
sd(unlist(sim1[1,3,]))/sqrt(n)
sd(unlist(sim1[1,4,]))/sqrt(n)
sd(unlist(sim1[1,5,]))/sqrt(n)

perm_sim1 <- sim1
save(perm_sim1, file = "simulations/perm_sim/permutation_res_plus.rda")


