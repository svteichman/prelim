library(e1071)

# Prep X and Y ----

## process HIV data from http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt to send to matlab

hivdata=read.table('genetic_application/NNRTI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)

y_cols=4:6;x_cols=7:246 # locations in the data matrix

is.letter <- function(x) grepl("[[:alpha:]]", x)

# remove the one sample with nonstandard mutation names
inc_err <- apply(hivdata[,x_cols], c(1:2), function(x) grepl("\\*", x)) 
row_err <- rowSums(inc_err) == 0
rmv <- which(row_err == FALSE)
hivdata=hivdata[-rmv,]

mutnums=read.table('genetic_application/NNRTI/treatment_mut_table.txt')

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

hivX <- X
Xnums <- Xnums
hivY <- Y

# Generate non-Gaussian errors ----
X <- hivX
Y <- hivY
freq <- which(colSums(X) > 2)
X <- X[,freq]

dupcols <- which(colSums(abs(cor(X)-1) < 1e-10) > 1)
if (length(dupcols > 0)) {
  X <- X[,-dupcols]
}

# make empirical distribution 
N <- nrow(X)
distr <- rep(0,3*N)
curr <- 1
for (i in 1:3) {
  rm <- which(Y[,i] == 999)
  y <- Y[-rm,i]
  x <- X[-rm,]
  
  dupcols <- which(colSums(abs(cor(x)-1) < 1e-10) > 1)
  if (length(dupcols > 0)) {
    x <- x[,-dupcols]
  }
  
  n <- nrow(x)
  proj <- diag(nrow = n) - x %*% solve(t(x) %*% x) %*% t(x) 
  vals <- proj %*% y
  distr[curr:(curr + n-1)] <- vals
  curr <- curr + n
}
distr <- distr[1:(curr-1)]
distr_sc <- scale(distr, center = T, scale = T)

save(X, file = "genetic_application/non_gaussian_sim/X.rda")
save(distr_sc, file = "genetic_application/non_gaussian_sim/distr.rda")


