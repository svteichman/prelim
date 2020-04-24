old_pi_X <- read.table("data/PI_X.txt")
old_pi_Y <- read.table("data/PI_Y.txt")
old_pi_Xnum <- read.table("data/PI_Xnums.txt")
old_nrti_X <- read.table("data/NRTI_X.txt")
old_nrti_Y <- read.table("data/NRTI_Y.txt")
old_nrti_Xnum <- read.table("data/NRTI_Xnums.txt")
old_nnrti_X <- read.table("data/NNRTI_X.txt")
old_nnrti_Y <- read.table("data/NNRTI_Y.txt")
old_nnrti_Xnum <- read.table("data/NNRTI_Xnums.txt")

old_pi_Y <- old_pi_Y[,1:7] #removing TPV, DRV 
old_nrti_Y <- old_nrti_Y[,c(1,2,3,4,6,8)] #removing DDC, FTC 
old_nnrti_Y <- old_nnrti_Y[,1:3] #removing ETR, RPV

tsm_pi <- read.delim('data/TMS_PI_frompaper.txt', fill = T, header = F)
tsm_pi <- tsm_pi[-1,]
pi_true <- tsm_pi$V2[tsm_pi$V2 >= 1]

tsm_nrti <- read.delim('data/TMS_NRTI_frompaper.txt', fill = T, header = F)
tsm_nrti <- tsm_nrti[-1,]
nrti_true <- tsm_nrti$V2[tsm_nrti$V2 >= 1]

tsm_nnrti <- read.delim('data/TMS_NNRTI_frompaper.txt', fill = T, header = F)
tsm_nnrti <- tsm_nnrti[-1,]
nnrti_true <- tsm_nnrti$V2[tsm_nnrti$V2 >= 1]

# go through analysis from tutorial 

knockoff_and_bhq <- function (X, y, q) {
  # Log-transform the drug resistance measurements.
  #y = log(y) # already done 
  
  # Remove patients with missing measurements.
  #missing = is.na(y)
  missing = y == 999 #recode based on coding from data 
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  # Run the knockoff filter.
  knock.gen = function(x) create.fixed(x, method='equi')
  result = knockoff.filter(X, y, fdr=fdr, knockoffs=knock.gen, statistic=stat.glmnet_lambdasmax)
  #knockoff_selected = names(result$selected)
  select <- unlist(result$selected)
  knockoff_selected = old_pi_Xnum$V1[select]
  
  # Run BHq.
  p = ncol(X)
  lm.fit = lm(y ~ X - 1) # no intercept
  p.values = coef(summary(lm.fit))[,4]
  cutoff = max(c(0, which(sort(p.values) <= fdr * (1:p) / p)))
  #bhq_selected = names(which(p.values <= fdr * cutoff / p))
  select_bh = which(p.values <= fdr * cutoff / p)
  bhq_selected = old_pi_Xnum$V1[select_bh]
  
  list(Knockoff = knockoff_selected, BHq = bhq_selected)
}

fdr = 0.20
results = lapply(old_pi_Y, function(y) knockoff_and_bhq(as.matrix(old_pi_X), y, fdr))

#get_position <- function(x)
 # sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

comparisons <- lapply(results, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(selected) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, pi_true))
    list(true_discoveries = discoveries - false_discoveries,
         false_discoveries = false_discoveries,
         fdp = false_discoveries / max(1, discoveries))
  })
})
names(comparisons) <- c("APV","ATV","IDV","LPV","NFV","FTV","SQV")
for (drug in names(comparisons)) {
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}

png("pi_barplots.png")
par(mfrow = c(3,3))
for (drug in names(comparisons)) {
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}
dev.off()
