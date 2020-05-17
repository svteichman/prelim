library(prelim.knockoffs)
library(knockoff)
library(mvtnorm)

# Knockoff 

normalize <- function(X) {
  X.scaled <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  X_scaled <- X.scaled[, ]
  return(X_scaled)
}

generate_dat <- function(n = 3000, p = 1000, k = 30) {
  mean <- rep(0,p)
  cov <- diag(x = 1, nrow = p)
  X <- rmvnorm(n, mean, cov)
  X <- normalize(X)
  z <- rnorm(n, 0, 1) 
  true_var <- sample(1:p, k)
  amp <- 3.5 *sample(c(-1,1), k, replace = TRUE)
  get_xb <- function(vec, true_var, amp) {
    return(sum(vec[true_var]*amp))
  }
  xb <- apply(X, 1, function(x) get_xb(x, true_var, amp))
  y <- xb + z
  return(list(X = X,y = y,true_var = true_var))
}

get_fdr <- function(n = 3000, p = 1000, k = 30) {
  dat <- generate_dat(n, p, k)
  X <- dat$X
  y <- dat$y
  true_var <- dat$true_var
  knock_plus_equi <- knockoff_filter(X, y, "equi", fdr = .2, plus = TRUE, randomize = TRUE)
  knock_equi <-  knockoff_filter(X, y, "equi", fdr = .2, plus = FALSE, randomize = TRUE)
  knock_plus_sdp <- knockoff_filter(X, y, "sdp", fdr = .2, plus = TRUE, randomize = TRUE)
  knock_sdp <- knockoff_filter(X, y, "sdp", fdr = .2, plus = FALSE, randomize = TRUE)
  bhq_base <- run_BHq(X, y, .2, "base")
  bhq_log <- run_BHq(X, y, .2, "log")
  bhq_white_noise <- run_BHq(X, y, .2, "white_noise")
  
  knock.gen.equi <- function(x) create.fixed(x, method = "equi", randomize = TRUE)
  knock.gen.sdp <- function(x) create.fixed(x, method = "sdp", randomize = TRUE)
  knock_plus_equi_R <- knockoff.filter(X, y, fdr = .2, knockoffs=knock.gen.equi, 
                                       statistic=stat.glmnet_lambdasmax, offset = 1)
  knock_equi_R <-  knockoff.filter(X, y, fdr = .2, knockoffs=knock.gen.equi, 
                                   statistic=stat.glmnet_lambdasmax, offset = 0)
  knock_plus_sdp_R <- knockoff.filter(X, y, fdr = .2, knockoffs=knock.gen.sdp, 
                                      statistic=stat.glmnet_lambdasmax, offset = 1)
  knock_sdp_R <- knockoff.filter(X, y, fdr = .2, knockoffs=knock.gen.sdp, 
                                 statistic=stat.glmnet_lambdasmax, offset = 0)
  
  res <- data.frame(tmp = 1)
  calc_fdr <- function(selected, true_var) {
    fdr <- ifelse(length(selected)==0, 0, 1 - sum(selected %in% true_var)/max(1,length(selected)))
    return(fdr)
  }
  res$knock_plus_equi_fdr <- calc_fdr(knock_plus_equi$selected, true_var)
  res$tmp <- NULL
  res$knock_equi_fdr <- calc_fdr(knock_equi$selected, true_var)
  res$knock_plus_sdp_fdr <- calc_fdr(knock_plus_sdp$selected, true_var)
  res$knock_sdp_fdr <- calc_fdr(knock_sdp$selected, true_var)
  res$bhq_base_fdr <- calc_fdr(bhq_base$selected, true_var)
  res$bhq_log_fdr <- calc_fdr(bhq_log$selected, true_var)
  res$bhq_white_noise_fdr <- calc_fdr(bhq_white_noise$selected, true_var)
  
  res$knock_plus_equi_fdr_R <- calc_fdr(knock_plus_equi_R$selected, true_var)
  res$knock_equi_fdr_R <- calc_fdr(knock_equi_R$selected, true_var)
  res$knock_plus_sdp_fdr_R <- calc_fdr(knock_plus_sdp_R$selected, true_var)
  res$knock_sdp_fdr_R <- calc_fdr(knock_sdp_R$selected, true_var)
  
  calc_power <- function(selected, true_var) {
    return(sum(selected %in% true_var)/length(true_var))
  }
  res$knock_plus_equi_pwr <- calc_power(knock_plus_equi$selected, true_var)
  res$knock_equi_pwr <- calc_power(knock_equi$selected, true_var)
  res$knock_plus_sdp_pwr <- calc_power(knock_plus_sdp$selected, true_var)
  res$knock_sdp_pwr <- calc_power(knock_sdp$selected, true_var)
  res$bhq_base_pwr <- calc_power(bhq_base$selected, true_var)
  res$bhq_log_pwr <- calc_power(bhq_log$selected, true_var)
  res$bhq_white_noise_pwr <- calc_power(bhq_white_noise$selected, true_var)
  
  res$knock_plus_equi_pwr_R <- calc_power(knock_plus_equi_R$selected, true_var)
  res$knock_equi_pwr_R <- calc_power(knock_equi_R$selected, true_var)
  res$knock_plus_sdp_pwr_R <- calc_power(knock_plus_sdp_R$selected, true_var)
  res$knock_sdp_pwr_R <- calc_power(knock_sdp_R$selected, true_var)
  
  return(res)
}

bhq_res <- get_fdr()  

args <- commandArgs(trailingOnly = TRUE)
save(bhq_res, file = paste0("bhq_res",args,".rda"))


