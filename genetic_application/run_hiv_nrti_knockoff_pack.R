library(prelim.knockoffs)
library(knockoff)
library(tidyverse)
library(gridExtra)
library(grid)
library(lattice)

# Prep X and Y ----

## process HIV data from http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt to send to matlab

hivdata=read.table('genetic_application/NRTI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)

y_cols=4:9;x_cols=10:249 # locations in the data matrix

is.letter <- function(x) grepl("[[:alpha:]]", x)

# remove the one sample with nonstandard mutation names
inc_err <- apply(hivdata[,x_cols], c(1:2), function(x) grepl("\\*", x)) 
row_err <- rowSums(inc_err) == 0
rmv <- which(row_err == FALSE)
hivdata=hivdata[-rmv,]

mutnums=read.table('genetic_application/NRTI/treatment_mut_table.txt')
nrti_tsm <- read.delim(file = "genetic_application/NRTI/NRTI_tsm.txt", header = FALSE)

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

# run analysis -----
run_hiv_data <- function(Xthresh,q) {
  ny <- dim(hivY)[2]
  
  tmtpos <- unique(nrti_tsm[,1])
  Xpos <- Xnums[,1]
  
  results <- array(0, dim = c(ny,3,2))
  for (iy in 1:ny) {
    dataloc <- which(hivY[,iy] != 999)
    Y <- hivY[dataloc,iy]
    X <- hivX[dataloc,]
    keepcols <- which(colSums(X) >= Xthresh)
    X <- X[,keepcols]
    
    #adding thing about correlation
    dupcols <- which(colSums(abs(cor(X)-1) < 1e-4) > 1)
    if (length(dupcols > 0)) {
      keepcols <- keepcols[-dupcols]
      X <- X[,-dupcols]
    }
    
    n <- nrow(X)
    p <- ncol(X)
    knock.gen.sdp <- function(x) create.fixed(x, method = "sdp", y = Y, randomize = TRUE)
    
    #S_knock <- keepcols[knockoff_filter(X,Y,method = "sdp",fdr= q, plus = FALSE, 
    #                                    randomize = TRUE)$selected]
    S_knock <- keepcols[knockoff.filter(X, Y, fdr = .2, knockoffs=knock.gen.sdp, 
                    statistic=stat.glmnet_lambdasmax, offset = 0)$selected]
    P_knock <- unique(Xpos[S_knock])
    #S_kplus <- keepcols[knockoff_filter(X,Y,method = "sdp",fdr= q, plus = TRUE, 
    #                                    randomize = TRUE)$selected]
    S_kplus <- keepcols[knockoff.filter(X, Y, fdr = .2, knockoffs=knock.gen.sdp, 
                    statistic=stat.glmnet_lambdasmax, offset = 1)$selected]
    P_kplus <- unique(Xpos[S_kplus])
    
    S_bhq <- keepcols[run_BHq(X,Y,q,"base")$selected]
    P_bhq <- unique(Xpos[S_bhq])
    
    results[iy,1,] <- c(length(P_knock), sum(P_knock %in% tmtpos))
    results[iy,2,] <- c(length(P_kplus), sum(P_kplus %in% tmtpos))
    results[iy,3,] <- c(length(P_bhq), sum(P_bhq %in% tmtpos))
    
  }
  return(results)
}

res_100 <- replicate(100, run_hiv_data(3,.2))

save(res_100, file = "genetic_application/nrti_res_100_knock_pack.rda")
