library(knockoff)

## process HIV data from http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt to send to matlab

hivdata=read.table('PI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)

y_cols=4:10;x_cols=11:109 # locations in the data matrix

is.letter <- function(x) grepl("[[:alpha:]]", x)

rmv=which(rowSums(hivdata[,x_cols]=='Krk')==1)
# remove the one sample with nonstandard mutation names
hivdata=hivdata[-rmv,]

#mutnums=read.table('~/Dropbox/matlab/HIVdata/mut_to_num.txt')
mutnums=read.table('PI/treatment_mut_table.txt')

X=Xnames=Xnums=NULL;n=dim(hivdata)[1]
for(i in x_cols){
  Xi=NULL
  col=hivdata[,i]
  mutnames=NULL
  for(j in 1:n){
    muts=strsplit(toString(col[j]),'')[[1]]
    for(k in 1:length(muts)){
      if(is.letter(muts[k])){
        if(any(mutnames==muts[k])){
          l=which(mutnames==muts[k])
          Xi[j,l]=1
        }else{
          vec=rep(0,n);vec[j]=1
          Xi=cbind(Xi,vec)
          mutnames=c(mutnames,muts[k])
          mutnum=mutnums[which(mutnums[,1]==muts[k]),2]
          Xnames=rbind(Xnames,c(as.numeric(strsplit(names(hivdata)[i],'P')[[1]][2]),muts[k]))
          Xnums=rbind(Xnums,c(as.numeric(strsplit(names(hivdata)[i],'P')[[1]][2]),mutnum))
        }
      }
    }
  }
  X=cbind(X,Xi)
}

Y=NULL
for(i in y_cols){
  find_na=which(is.na(hivdata[,i]))
  col=log(hivdata[,i]);col[find_na]=999
  Y=cbind(Y,col)
}

#### DONE WITH SECTION FROM RINA'S CODE - USE CODE TO REPRODUCE RESULTS ! 

pi_tsm <- read.delim(file = "data/PI_tsm.txt", header = FALSE, stringsAsFactors = FALSE)
names(pi_tsm) = c('Position', 'Mutations')
nrti_tsm <- read.delim(file = "data/NRTI_tsm.txt", header = FALSE, stringsAsFactors = FALSE)
names(nrti_tsm) = c('Position', 'Mutations')
nnrti_tsm <- read.delim(file = "data/NNRTI_tsm.txt", header = FALSE, stringsAsFactors = FALSE)
names(nnrti_tsm) = c('Position', 'Mutations')

tsm_pi <- read.delim('data/TMS_PI_frompaper.txt', fill = T, header = F)
tsm_pi <- tsm_pi[-1,]
pi_true <- tsm_pi$V2[tsm_pi$V2 >= 1]
tsm_nrti <- read.delim('data/TMS_NRTI_frompaper.txt', fill = T, header = F)
tsm_nrti <- tsm_nrti[-1,]
nrti_true <- tsm_nrti$V2[tsm_nrti$V2 >= 1]
tsm_nnrti <- read.delim('data/TMS_NNRTI_frompaper.txt', fill = T, header = F)
tsm_nnrti <- tsm_nnrti[-1,]
nnrti_true <- tsm_nnrti$V2[tsm_nnrti$V2 >= 1]

knockoff_and_bhq <- function (X, y, q, Xnums) {
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
  knock.gen = function(x) create.fixed(x, method='equi', randomize = F)
  result = knockoff.filter(X, y, fdr=fdr, knockoffs=knock.gen, statistic=stat.glmnet_lambdasmax)
  #knockoff_selected = names(result$selected)
  select <- unlist(result$selected)
  knockoff_selected = Xnums[select]
  
  # Run BHq.
  p = ncol(X)
  lm.fit = lm(y ~ X - 1) # no intercept
  p.values = coef(summary(lm.fit))[,4]
  cutoff = max(c(0, which(sort(p.values) <= fdr * (1:p) / p)))
  #bhq_selected = names(which(p.values <= fdr * cutoff / p))
  select_bh = which(p.values <= fdr * cutoff / p)
  bhq_selected = Xnums[select_bh]
  
  list(Knockoff = knockoff_selected, BHq = bhq_selected)
}
colnames(Y) <- c("APV","ATV","IDV","LPV","NFV","RTV","SQV")
fdr = 0.20
results = apply(Y, 2, function(y) knockoff_and_bhq(X, y, fdr, Xnums))

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

names(comparisons) <- c("APV","ATV","IDV","LPV","NFV","RTV","SQV")
for (drug in names(comparisons)) {
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}


png("pi_barplots_2.png")
par(mfrow = c(3,3))
for (drug in names(comparisons)) {
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}
dev.off()


# RUN 10 TIMES!!!!
knock_true <- array(data = 0, dim = c(100,7))
knock_false <- array(data = 0, dim = c(100,7))
bhq_true <- array(data = 0, dim = c(100,7))
bhq_false <- array(data = 0, dim = c(100,7))
#results = apply(Y, 2, function(y) knockoff_and_bhq(X, y, fdr, Xnums))
#names(comparisons) <- c("APV","ATV","IDV","LPV","NFV","RTV","SQV")
for (i in 1:100) {
  results = apply(Y, 2, function(y) knockoff_and_bhq(X, y, fdr, Xnums))
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
  #print(results$APV[1:5])
  for (j in 1:7) {
    drug <- colnames(Y)[j]
    plot_data = do.call(cbind, comparisons[[drug]])
    plot_data = plot_data[c('true_discoveries','false_discoveries'),]
    knock_true[i,j] <- unlist(plot_data[1,1])
    knock_false[i,j] <- unlist(plot_data[2,1])
    bhq_true[i,j] <- unlist(plot_data[1,2])
    bhq_false[i,j] <- unlist(plot_data[2,2])
  }
}

knock_true_avg <- colMeans(knock_true)
knock_false_avg <- colMeans(knock_false)
bhq_true_avg <- colMeans(bhq_true)
bhq_false_avg <- colMeans(bhq_false)

png("pi_barplots_100_mean.png")
par(mfrow = c(3,3))
for (j in 1:7) {
  drug <- names(comparisons)[j]
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  plot_data[1,1] <- knock_true_avg[j]
  plot_data[1,2] <- bhq_true_avg[j]
  plot_data[2,1] <- knock_false_avg[j]
  plot_data[2,2] <- bhq_false_avg[j]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}
dev.off()

knock_true_quant <- apply(knock_true, 2, function(p) quantile(p))
knock_false_quant <- apply(knock_false, 2, function(p) quantile(p))
bhq_true_quant <- apply(bhq_true, 2, function(p) quantile(p))
bhq_false_quant <- apply(bhq_false, 2, function(p) quantile(p))

png("pi_barplots_n100_p25.png")
par(mfrow = c(3,3))
for (j in 1:7) {
  drug <- names(comparisons)[j]
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  plot_data[1,1] <- knock_true_quant[2,j]
  plot_data[1,2] <- bhq_true_quant[2,j]
  plot_data[2,1] <- knock_false_quant[2,j]
  plot_data[2,2] <- bhq_false_quant[2,j]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}
dev.off()

png("pi_barplots_n100_p50.png")
par(mfrow = c(3,3))
for (j in 1:7) {
  drug <- names(comparisons)[j]
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  plot_data[1,1] <- knock_true_quant[3,j]
  plot_data[1,2] <- bhq_true_quant[3,j]
  plot_data[2,1] <- knock_false_quant[3,j]
  plot_data[2,2] <- bhq_false_quant[3,j]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}
dev.off()

png("pi_barplots_n100_p75.png")
par(mfrow = c(3,3))
for (j in 1:7) {
  drug <- names(comparisons)[j]
  plot_data = do.call(cbind, comparisons[[drug]])
  plot_data = plot_data[c('true_discoveries','false_discoveries'),]
  plot_data[1,1] <- knock_true_quant[4,j]
  plot_data[1,2] <- bhq_true_quant[4,j]
  plot_data[2,1] <- knock_false_quant[4,j]
  plot_data[2,2] <- bhq_false_quant[4,j]
  barplot(as.matrix(plot_data), main = paste('Resistance to', drug),
          col = c('navy','orange'), ylim = c(0,40))
}
dev.off()