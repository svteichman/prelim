library(tidyverse)
library(gridExtra)
library(grid)
library(lattice)

#PI

load(file = "genetic_application/rda_files/pi_res_100.rda")
pi_res <- res_100

load(file = "genetic_application/rda_files/pi_res_100_knock_pack.rda")
pi_res_ko <- res_100

hivdata=read.table('genetic_application/PI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)
y_cols=4:10;x_cols=11:109 # locations in the data matrix
pi_tsm <- read.delim(file = "genetic_application/PI/PI_tsm.txt", header = FALSE)
load("genetic_application/PI/hivX.rda")
load("genetic_application/PI/hivY.rda")

#info to add to plots 
drugnames <- names(hivdata[y_cols])
num_drug <- length(drugnames)

tot_tsm <- length(unique(pi_tsm$V1))

size_res <- matrix(nrow = num_drug, ncol = 2)
ny <- num_drug 
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

plot_fun <- function(tmp1, tmp2, file_name) {
  file = paste0("genetic_application/figs/", file_name)
  plot_res <- data.frame(method = rep(c(rep("Knockoff",num_drug),rep("Knockoff (R)",num_drug)),2),
                         value = c(tmp1[,1,2],tmp2[,1,2],tmp1[,1,1]-tmp1[,1,2],tmp2[,1,1]-tmp2[,1,2]),
                         type = c(rep("In TSM list",2*num_drug), rep("Not in TSM list",2*num_drug)),
                         drug = rep(drugnames,4)) 
  plot_res$type <- factor(plot_res$type, levels = c("In TSM list","Not in TSM list"))
  plot_res$method <- factor(plot_res$method, levels = c("Knockoff","Knockoff (R)"))
  labs <- paste0("Resistance to ",drugnames)
  plot_res$drug <- factor(plot_res$drug, labels = labs)
  
  dat_text <- data.frame(label = footers, drug = drugnames, x = rep(0.5,7),
                         y = rep(0.9,7))
  
  ggplot(plot_res, aes(x = method, y = value,  fill = type)) + 
    geom_col(position = position_stack(reverse = TRUE), width = .7) +
    facet_wrap(~drug, scales = "free") + ylim(c(0,40)) + 
    scale_fill_manual(values = c("midnightblue","chocolate2")) +
    ylab("# HIV-1 protease positions selected") + xlab(element_blank()) + 
    geom_hline(yintercept = tot_tsm, linetype = "dashed") +
    theme_bw()  +
    theme(legend.title = element_blank(), legend.position = c(0.5, 0.2),
          panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_text(data=data.frame(x= "Knockoff", y=37, label=footers, 
                              drug=labs), 
              aes(x = x, y = y, label=label), inherit.aes=FALSE, size = 3, nudge_x = .5 )
  ggsave(file = file)
}

# mean
res_avg1 <- apply(pi_res, 1:3, mean)
res_avg2 <- apply(pi_res_ko, 1:3, mean)
plot_fun(res_avg1, res_avg2, "pi_avg_compare.png")


#NRTI

load(file = "genetic_application/rda_files/nrti_res_100.rda")
nrti_res <- res_100

load(file = "genetic_application/rda_files/nrti_res_100_knock_pack.rda")
nrti_res_ko <- res_100

hivdata=read.table('genetic_application/NRTI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)
y_cols=4:9;x_cols=10:249 # locations in the data matrix
nrti_tsm <- read.delim(file = "genetic_application/NRTI/NRTI_tsm.txt", header = FALSE)
load("genetic_application/NRTI/hivX.rda")
load("genetic_application/NRTI/hivY.rda")

#info to add to plots 
drugnames <- names(hivdata[y_cols])
num_drug <- length(drugnames)

tot_tsm <- length(unique(nrti_tsm$V1))

size_res <- matrix(nrow = num_drug, ncol = 2)
ny <- num_drug 
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

plot_fun <- function(tmp1, tmp2, file_name) {
  file = paste0("genetic_application/figs/", file_name)
  plot_res <- data.frame(method = rep(c(rep("Knockoff",num_drug),rep("Knockoff (R)",num_drug)),2),
                         value = c(tmp1[,1,2],tmp2[,1,2],tmp1[,1,1]-tmp1[,1,2],tmp2[,1,1]-tmp2[,1,2]),
                         type = c(rep("In TSM list",2*num_drug), rep("Not in TSM list",2*num_drug)),
                         drug = rep(drugnames,4)) 
  plot_res$type <- factor(plot_res$type, levels = c("In TSM list","Not in TSM list"))
  plot_res$method <- factor(plot_res$method, levels = c("Knockoff","Knockoff (R)"))
  labs <- paste0("Resistance to ",drugnames)
  plot_res$drug <- factor(plot_res$drug, labels = labs)
  
  dat_text <- data.frame(label = footers, drug = drugnames, x = rep(0.5,num_drug),
                         y = rep(0.9,num_drug))
  
  ggplot(plot_res, aes(x = method, y = value,  fill = type)) + 
    geom_col(position = position_stack(reverse = TRUE), width = .7) +
    facet_wrap(~drug, scales = "free") + ylim(c(0,31)) + 
    scale_fill_manual(values = c("midnightblue","chocolate2")) +
    ylab("# HIV-1 RT positions selected") + xlab(element_blank()) + 
    geom_hline(yintercept = tot_tsm, linetype = "dashed") +
    theme_bw()  +
    theme(legend.title = element_blank(), legend.position = "none",
          panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_text(data=data.frame(x= "Knockoff", y=27, label=footers, 
                              drug=labs), 
              aes(x = x, y = y, label=label), inherit.aes=FALSE, size = 3, nudge_x = .5 )
  ggsave(file = file, width = 6.4, height = 3.26)
}

# mean
res_avg1 <- apply(nrti_res, 1:3, mean)
res_avg2 <- apply(nrti_res_ko, 1:3, mean)
plot_fun(res_avg1, res_avg2, "nrti_avg_compare.png")




#NNRTI

load(file = "genetic_application/rda_files/nnrti_res_100.rda")
nnrti_res <- res_100

load(file = "genetic_application/rda_files/nnrti_res_100_knock_pack.rda")
nnrti_res_ko <- res_100

hivdata=read.table('genetic_application/NNRTI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)
y_cols=4:6;x_cols=7:246 # locations in the data matrix
nnrti_tsm <- read.delim(file = "genetic_application/NNRTI/NNRTI_tsm.txt", header = FALSE)
load("genetic_application/NNRTI/hivX.rda")
load("genetic_application/NNRTI/hivY.rda")

#info to add to plots 
drugnames <- names(hivdata[y_cols])
num_drug <- length(drugnames)

tot_tsm <- length(unique(nnrti_tsm$V1))

size_res <- matrix(nrow = num_drug, ncol = 2)
ny <- num_drug 
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

plot_fun <- function(tmp1, tmp2, file_name) {
  file = paste0("genetic_application/figs/", file_name)
  plot_res <- data.frame(method = rep(c(rep("Knockoff",num_drug),rep("Knockoff (R)",num_drug)),2),
                         value = c(tmp1[,1,2],tmp2[,1,2],tmp1[,1,1]-tmp1[,1,2],tmp2[,1,1]-tmp2[,1,2]),
                         type = c(rep("In TSM list",2*num_drug), rep("Not in TSM list",2*num_drug)),
                         drug = rep(drugnames,4)) 
  plot_res$type <- factor(plot_res$type, levels = c("In TSM list","Not in TSM list"))
  plot_res$method <- factor(plot_res$method, levels = c("Knockoff","Knockoff (R)"))
  labs <- paste0("Resistance to ",drugnames)
  plot_res$drug <- factor(plot_res$drug, labels = labs)
  
  dat_text <- data.frame(label = footers, drug = drugnames, x = rep(0.5,num_drug),
                         y = rep(0.9,num_drug))
  
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
res_avg1 <- apply(nnrti_res, 1:3, mean)
res_avg2 <- apply(nnrti_res_ko, 1:3, mean)
plot_fun(res_avg1, res_avg2, "nnrti_avg_compare.png")
