library(prelim.knockoffs)
library(mvtnorm)

# Knockoff 

normalize <- function(X) {
  X.scaled <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  X_scaled <- X.scaled[, ]
  return(X_scaled)
}

generate_dat <- function(n = 3000, p = 1000, k = 30, rho = 0) {
  mean <- rep(0,p)
  cov <- diag(x = 1, nrow = p)
  for (i in 1:p) {
    for (j in 1:p) {
      cov[i,j] <- rho^(abs(j-i))
    }
  }
  X <- rmvnorm(n, mean, cov)
  X <- normalize(X)
  z <- rnorm(n, 0, 1) 
  true_var <- sample(1:p, k)
  A <- 3.5 
  amp <- A *sample(c(-1,1), k, replace = TRUE)
  get_xb <- function(vec, true_var, amp) {
    return(sum(vec[true_var]*amp))
  }
  xb <- apply(X, 1, function(x) get_xb(x, true_var, amp))
  y <- xb + z
  return(list(X = X,y = y,true_var = true_var))
}

get_fdr <- function(n = 3000, p = 1000, k = 30, rho = 0) {
  dat <- generate_dat(n, p, k, rho)
  X <- dat$X
  y <- dat$y
  true_var <- dat$true_var
  knock_plus_sdp <- knockoff_filter(X, y, "sdp", fdr = .2, plus = TRUE)
  knock_sdp <- knockoff_filter_from_stats(X, y, fdr = .2, plus = FALSE, Xk = knock_plus_sdp$Xk,
                                          W = knock_plus_sdp$statistic)
  bhq_base <- run_BHq(X, y, .2, "base")
  
  res <- data.frame(tmp = 1)
  calc_fdr <- function(selected, true_var) {
    fdr <- ifelse(length(selected)==0, 0, 1 - sum(selected %in% true_var)/max(1,length(selected)))
    return(fdr)
  }
  res$knock_plus_sdp_fdr <- calc_fdr(knock_plus_sdp$selected, true_var)
  res$tmp <- NULL
  res$knock_sdp_fdr <- calc_fdr(knock_sdp$selected, true_var)
  res$bhq_base_fdr <- calc_fdr(bhq_base$selected, true_var)
  
  calc_power <- function(selected, true_var) {
    return(sum(selected %in% true_var)/length(true_var))
  }
  res$knock_plus_sdp_pwr <- calc_power(knock_plus_sdp$selected, true_var)
  res$knock_sdp_pwr <- calc_power(knock_sdp$selected, true_var)
  res$bhq_base_pwr <- calc_power(bhq_base$selected, true_var)
  
  return(res)
}

rho_seq <- seq(0, 0.9, by = .1)
bhq_sim2c <- get_fdr(rho = rho_seq[1])

for (i in 2:length(rho_seq)) {
  new_res <- get_fdr(rho = rho_seq[i])
  bhq_sim2c <- rbind(bhq_sim2c, new_res)
}

args <- commandArgs(trailingOnly = TRUE)
save(bhq_sim2c, file = paste0("BHq_sim2c",args,".rda"))

