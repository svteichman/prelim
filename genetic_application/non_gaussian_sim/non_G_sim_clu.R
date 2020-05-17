library(prelim.knockoffs)

generate_y <- function(X, distr) {
  N <- nrow(X)
  p <- ncol(X)
  true_var <- sample(1:p, 20)
  beta_j <- 3.5*rnorm(20, 0, 1)
  beta <- rep(0,p)
  beta[true_var] <- beta_j
  
  ind <- sample(1:length(distr), N, replace = TRUE)
  y <- X %*% beta + distr[ind]
  return(list(X = X, y = y, true_var = true_var))
}

get_fdr <- function(X, distr_sc) {
  dat <- generate_y(X, distr_sc)
  X <- dat$X
  
  dupcols <- which(colSums(abs(cor(X)-1) < 1e-4) > 1)
  if (length(dupcols > 0)) {
    X <- X[,-dupcols]
  }
  
  y <- dat$y
  true_var <- dat$true_var
  
  knock_plus_equi <- knockoff_filter(X, y, "equi", fdr = .2, plus = TRUE)
  knock_equi <- knockoff_filter(X, y, "equi", fdr = .2, plus = FALSE)
  knock_plus_sdp <- knockoff_filter(X, y, "sdp", fdr = .2, plus = TRUE)
  knock_sdp <- knockoff_filter(X, y, "sdp", fdr = .2, plus = FALSE)
  bhq_base <- run_BHq(X, y, fdr = .2, method = "base")
  
  res <- data.frame(tmp = 1)
  calc_fdr <- function(selected, true_var) {
    fdr <- ifelse(length(selected)==0, 0, 1 - sum(selected %in% true_var)/max(1,length(selected)))
    return(fdr)
  }
  res$knock_plus_equ_fdr <- calc_fdr(knock_plus_equi$selected, true_var)
  res$tmp <- NULL
  res$knock_equ_fdr <- calc_fdr(knock_equi$selected, true_var)
  res$knock_plus_sdp_fdr <- calc_fdr(knock_plus_sdp$selected, true_var)
  res$knock_sdp_fdr <- calc_fdr(knock_sdp$selected, true_var)
  res$bhq_base_fdr <- calc_fdr(bhq_base$selected, true_var)
  
  calc_power <- function(selected, true_var) {
    return(sum(selected %in% true_var)/length(true_var))
  }
  res$knock_plus_equi_pwr <- calc_power(knock_plus_equi$selected, true_var)
  res$knock_equi_pwr <- calc_power(knock_equi$selected, true_var)
  res$knock_plus_sdp_pwr <- calc_power(knock_plus_sdp$selected, true_var)
  res$knock_sdp_pwr <- calc_power(knock_sdp$selected, true_var)
  res$bhq_base_pwr <- calc_power(bhq_base$selected, true_var)
  
  return(res)
}

load(file = "X.rda")
load(file = "distr.rda")

sim_res <- get_fdr(X, distr_sc)  

args <- commandArgs(trailingOnly = TRUE)
save(sim_res, file = paste0("non_G_sim",args,".rda"))