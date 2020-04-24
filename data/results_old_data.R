library(tidyverse)
library(readr)
library(knockoff)
setwd('~/git/prelim/data')

X_pi <- read.table("PI_X.txt")
Y_pi <- read.table("PI_Y.txt")
Xnums_pi <- read.table("PI_Xnums.txt")

tsm_pi <- read.delim("TMS_PI_frompaper.txt", header=F)
tsm_pi <- tsm_pi[-1,]
pi_true_pos <- tsm_pi$V2[tsm_pi$V2>=1]
tsm_nrti <- read.delim("TMS_NRTI_frompaper.txt", header=F)
tsm_nrti <- tsm_nrti[-1,]
nrti_true_pos <- tsm_nrti$V2[tsm_nrti$V2>=1]
tsm_nnrti <- read.delim("TMS_NNRTI_frompaper.txt", header=F)
tsm_nnrti <- tsm_nnrti[-1,]
nnrti_true_pos <- tsm_nnrti$V2[tsm_nnrti$V2>=1]

# older data 
#pr_dat <- read.delim(file = "data/PR_data.txt", na.string = c('NA', ''), stringsAsFactors = FALSE)
#rt_dat <- read.delim(file = "data/RT_data.txt", na.string = c('NA', ''), stringsAsFactors = FALSE)
#pi_tsm_old <- read.delim(file = "data/TMS_PI_frompaper.txt", header = FALSE, stringsAsFactors = FALSE)
#names(pi_tsm_old) = c('Position', 'Mutations')
#nrti_tsm_old <- read.delim(file = "data/TMS_NRTI_frompaper.txt", header = FALSE, stringsAsFactors = FALSE)
#nnrti_tsm_old <- read.delim(file = "data/TMS_NNRTI_frompaper.txt", header = FALSE, stringsAsFactors = FALSE)

# read in data for 3 different types of drugs 
pi <- read.delim(file = "PI_DATA.txt", na.string = c('NA', ''), stringsAsFactors = FALSE)
nrti <- read.delim(file = "NRTI_DATA.txt", na.string = c('NA', ''), stringsAsFactors = FALSE)
nnrti <- read.delim(file = "NNRTI_DATA.txt", na.string = c('NA', ''), stringsAsFactors = FALSE)

# read in tsm files for each drug 
#pi_tsm <- read.delim(file = "data/PI_tsm.txt", header = FALSE, stringsAsFactors = FALSE)
#names(pi_tsm) = c('Position', 'Mutations')
#nrti_tsm <- read.delim(file = "data/NRTI_tsm.txt", header = FALSE, stringsAsFactors = FALSE)
#names(nrti_tsm) = c('Position', 'Mutations')
#nnrti_tsm <- read.delim(file = "data/NNRTI_tsm.txt", header = FALSE, stringsAsFactors = FALSE)
#names(nnrti_tsm) = c('Position', 'Mutations')

# function to remove rows with errors and non-standard mutations, * or lowercase letters (besides i,d)
rows_to_keep <- function(dat) {
  m_start <- which(names(dat) == "P1")
  m_end <- ncol(dat) 
  inc_err <- apply(dat[m_start:m_end], c(1:2), function(x) grepl("\\*|[a-d]|[e-h]|[j-z]", x)) 
  row_err <- rowSums(inc_err) == 0
}

# function to go from data frame to design matrix 
# variable for every combination of mutation {A,...,Z,i,d} and position 
get_design_mat <- function(dat) {
  n <- nrow(dat)
  m_start <- which(names(dat) == "P1")
  m_end <- ncol(dat) 
  mut <- c(LETTERS,"i","d")
  p <- length(mut)*(m_end-m_start+1)
  col_name_comp <- expand.grid(1:(m_end-m_start+1),mut,stringsAsFactors = F)
  col_names <- paste0("P",col_name_comp$Var1,".",col_name_comp$Var2)
  X <- matrix(0, n, p)
  for (i in 1:n) {
    mut_pos <- which(grepl("-",dat[i,m_start:m_end])==FALSE)
    mut_num <- length(mut_pos)
    if (mut_num > 0) {
      for (j in 1:mut_num) {
        mut_lett <- dat[i,mut_pos[j]-1+m_start]
        if (nchar(mut_lett) > 1) {
          mut_lett <- unlist(strsplit(dat[i,mut_pos[j]-1+m_start],""))
        }
        for (l in 1:length(mut_lett)) {
          col <- which(col_names == paste0("P",mut_pos[j],".",mut_lett[l]))
          X[i,col] <- 1 
        }
      }
    }
  }
  colnames(X) <- col_names
  return(X)
}

# wrapper function to go from data frame to X and Y 
get_dat <- function(dat) {
  valid_rows <- rows_to_keep(dat)
  dat <- dat[valid_rows,]
  X <- get_design_mat(dat) 
  mode(X) <- 'numeric'
  X <- X[,colSums(X) > 0]
  m_start <- which(names(dat) == "P1")
  Y <- dat[,4:(m_start-1)]
  return(list(X,Y))
}

# do knockoffs and BH procedure on data 
select_var <- function(X,y, fdr) {
  y <- log(y) 
  # remove missing
  y_miss <- is.na(y)
  y <- y[!y_miss]
  X <- X[!y_miss,]
  # remove columns with few (< 3) observations
  X_sparse <- colSums(X) < 3
  X <- X[,!X_sparse]
  # remove approximately duplicated columns 
  X_dup <- colSums(abs(cor(X)) > (1-1e-4)) > 1
  X <- X[,!X_dup] 
  
  # run knockoffs (fixed-x, equicorrelated method of choosing s)
  knock = function(x) create.fixed(x, method='equi', y = y)
  result = knockoff.filter(X, y, fdr=fdr, knockoffs=knock, statistic=stat.glmnet_lambdasmax)
  knockoff_selected = names(result$selected)
  
  # Run BHq.
  p = ncol(X)
  mod = lm(y ~ X - 1) # remove intercept 
  p.values = coef(summary(mod))[,4]
  cutoff = max(c(0, which(sort(p.values) <= fdr * (1:p) / p)))
  bhq_selected = substring(names(which(p.values <= fdr * cutoff / p)),2)
  
  return(list(Knockoff = knockoff_selected, BHq = bhq_selected))
}

# function to get numbers of true/false discoveries 
get_res <- function(dat, tsm) {
  x_y <- get_dat(dat)
  X <- x_y[[1]]
  Y <- x_y[[2]]
  fdr <- 0.20
  results <- lapply(Y, function(y) select_var(X, y, fdr))
  drug_num <- length(results)
  for (i in 1:drug_num) {
    results[[i]] <- lapply(results[[i]], function(res) sort(unique(parse_number(res))))
  }
  true <- tsm
  res <- lapply(results, function(drug) {
    lapply(drug, function(vars) {
      disc <- length(vars)
      false <- length(setdiff(vars, true))
      list(true_discoveries = disc - false,
           false_discoveries = false,
           fdr = false/max(1,disc))
    })
  })
}

# Plot results for PI type drugs 
#pi_res <- get_res(pi, pi_tsm)
#true_num <- length(pi_tsm$Position)
pi_res <- get_res(pi, pi_true_pos)
true_num <- length(pi_true_pos)
drugs <- names(pi_res)
for (drug in drugs) {
  plot_df <- matrix(unlist(pi_res[[drug]]),ncol = 2)[-3,]
  colnames(plot_df) <- c("Knockoff","BHq")
  barplot(plot_df, main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
  abline(h = true_num, lty = 2)
}

# Plot results for NRTI type drugs 
#nrti_res <- get_res(nrti, nrti_tsm)
#true_num <- length(nrti_tsm$Position)
nrti_res <- get_res(nrti, nrti_true_pos)
true_num <- length(nrti_true_pos)
drugs <- names(nrti_res)
for (drug in drugs) {
  plot_df <- matrix(unlist(nrti_res[[drug]]),ncol = 2)[-3,]
  colnames(plot_df) <- c("Knockoff","BHq")
  barplot(plot_df, main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
  abline(h = true_num, lty = 2)
}

# Plot results for NNRTI type drugs 
#nnrti_res <- get_res(nnrti, nnrti_tsm)
#true_num <- length(nnrti_tsm$Position)
nnrti_res <- get_res(nnrti, nnrti_true_pos)
true_num <- length(nnrti_true_post)
drugs <- names(nnrti_res)
for (drug in drugs) {
  plot_df <- matrix(unlist(nnrti_res[[drug]]),ncol = 2)[-3,]
  colnames(plot_df) <- c("Knockoff","BHq")
  barplot(plot_df, main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
  abline(h = true_num, lty = 2)
}
