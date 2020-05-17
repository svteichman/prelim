library(prelim.knockoffs)
library(tidyverse)
library(gridExtra)
library(grid)
library(lattice)

# Prep X and Y ----
# Part of Rina's code to process data, received over email 
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
nnrti_tsm <- read.delim(file = "genetic_application/NNRTI/NNRTI_tsm.txt", header = FALSE)

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

save(hivX, file = "genetic_application/NNRTI/hivX.rda")
save(hivY, file = "genetic_application/NNRTI/hivY.rda")
save(Xnums, file = "genetic_application/NNRTI/Xnums.rda")

# run analysis -----
# my code to run analysis
run_hiv_data <- function(Xthresh,q) {
  ny <- dim(hivY)[2]
  
  tmtpos <- unique(nnrti_tsm[,1])
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
    S_knock <- keepcols[knockoff_filter(X,Y,method = "sdp",fdr= q, plus = FALSE, 
                                        randomize = TRUE)$selected]
    P_knock <- unique(Xpos[S_knock])
    S_kplus <- keepcols[knockoff_filter(X,Y,method = "sdp",fdr= q, plus = TRUE, 
                                        randomize = TRUE)$selected]
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

save(res_100, file = "genetic_application/rda_files/nnrti_res_100.rda")
load(file = "genetic_application/rda_files/nnrti_res_100.rda")

#info to add to plots 
drugnames <- names(hivdata[4:6])
num_drug <- length(drugnames)

tot_tsm <- length(unique(nnrti_tsm$V1))

size_res <- matrix(nrow = num_drug, ncol = 2)
ny <- dim(hivY)[2]
for (iy in 1:ny) {
  dataloc <- which(hivY[,iy] != 999)
  Y <- hivY[dataloc,iy]
  X <- hivX[dataloc,]
  keepcols <- which(colSums(X) >= 3)
  X <- X[,keepcols]
  
  #adding thing about correlation
  dupcols <- which(colSums(abs(cor(X)-1) < 1e-4) > 1)
  if (length(dupcols > 0)) {
    keepcols <- keepcols[-dupcols]
    X <- X[,-dupcols]
  }
  
  size_res[iy,1] <- nrow(X)
  size_res[iy,2] <- ncol(X)
}
footers <- paste0("n=",size_res[,1],", p=",size_res[,2])

# plot results 

plot_fun <- function(tmp, file_name) {
  file = paste0("genetic_application/figs/", file_name)
  plot_res <- data.frame(method = rep(c(rep("Knockoff",num_drug),rep("Knockoff+",num_drug), rep("BHq",num_drug)),2),
                         value = c(tmp[,1,2],tmp[,2,2],tmp[,3,2],tmp[,1,1]-tmp[,1,2],tmp[,2,1]-tmp[,2,2],tmp[,3,1]-tmp[,3,2]),
                         type = c(rep("In TSM list",3*num_drug), rep("Not in TSM list",3*num_drug)),
                         drug = rep(drugnames,6)) %>%
    filter(method != "Knockoff+") 
  plot_res$type <- factor(plot_res$type, levels = c("In TSM list","Not in TSM list"))
  plot_res$method <- factor(plot_res$method, levels = c("Knockoff","BHq"))
  labs <- paste0("Resistance to ",drugnames)
  plot_res$drug <- factor(plot_res$drug, levels = drugnames, labels = labs)
  
  ggplot(plot_res, aes(x = method, y = value,  fill = type)) + 
    geom_col(position = position_stack(reverse = TRUE), width = .7) +
    facet_wrap(~drug, scales = "free") + ylim(c(0,35)) + 
    scale_fill_manual(values = c("midnightblue","chocolate2")) +
    ylab("# positions selected") + xlab(element_blank()) + 
    geom_hline(yintercept = tot_tsm, linetype = "dashed") +
    theme_bw()  +
    theme(legend.title = element_blank(), legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_text(data=data.frame(x= "Knockoff", y=32, label=footers, 
                              drug=labs), 
              aes(x = x, y = y, label=label), inherit.aes=FALSE, size = 3, nudge_x = .5 )
  ggsave(file = file, width = 6.4, height = 1.63)
}

# mean 
avg_res <- apply(res_100, 1:3, mean)
plot_fun(avg_res, "nnrti_avg.png")


plot_quant <- function(tmp1, tmp2, file_name) {
  file = paste0("genetic_application/figs/", file_name)
  plot_res <- data.frame(method = rep(c(rep("10th",num_drug),rep("90th",num_drug)),2),
                         value = c(tmp1[,1,2],tmp2[,1,2],tmp1[,1,1]-tmp1[,1,2],tmp2[,1,1]-tmp2[,1,2]),
                         type = c(rep("In TSM list",2*num_drug), rep("Not in TSM list",2*num_drug)),
                         drug = rep(drugnames,4)) 
  plot_res$type <- factor(plot_res$type, levels = c("In TSM list","Not in TSM list"))
  plot_res$method <- factor(plot_res$method, levels = c("10th","90th"))
  labs <- paste0("Resistance to ",drugnames)
  plot_res$drug <- factor(plot_res$drug, labels = labs)
  
  ggplot(plot_res, aes(x = method, y = value,  fill = type)) + 
    geom_col(position = position_stack(reverse = TRUE), width = .7) +
    facet_wrap(~drug, scales = "free") + ylim(c(0,40)) + 
    scale_fill_manual(values = c("midnightblue","chocolate2")) +
    ylab("# positions selected") + xlab("Quantiles of Knockoff Results") + 
    geom_hline(yintercept = tot_tsm, linetype = "dashed") +
    theme_bw()  +
    theme(legend.title = element_blank(), legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_text(data=data.frame(x= "10th", y=37, label=footers, 
                              drug=labs), 
              aes(x = x, y = y, label=label), inherit.aes=FALSE, size = 3, nudge_x = .5 )
  ggsave(file = file, width = 6.4, height = 1.63)
}

# Plotting 10th and 90th quantiles  
res_10_quant <- apply(res_100, 1:3, function(x) quantile(x, probs = .1))
res_90_quant <- apply(res_100, 1:3, function(x) quantile(x, probs = .9))
plot_quant(res_10_quant, res_90_quant, "nnrti_q10q90.png")


