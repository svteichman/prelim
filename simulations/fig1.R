library(prelim.knockoffs)
library(mvtnorm)
library(tidyverse)
library(latex2exp)

normalize <- function(X) {
  X.scaled <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  X_scaled <- X.scaled[, ]
  return(X_scaled)
}

generate_dat <- function() {
  n <- 200
  p <- 100
  mean <- rep(0,p)
  cov <- diag(x = 1, nrow = p)
  X <- rmvnorm(n, mean, cov)
  z <- rnorm(n, 0, .7) 
  x_20_sum <- apply(X, 1, function(x) sum(x[1:20]))
  y <- .14*x_20_sum + z
  return(list(X=X,y=y))
}

dat <- generate_dat()
X <- dat$X
y <- dat$y

res <- create_knockoffs(X, y, "equi")
knock <- res$Xk
X <- res$X
X_aug <- cbind(X, knock)

n <- nrow(X_aug)
p <- ncol(X_aug)

num_lambda <- 500
lambda_max <- max(abs(t(X_aug) %*% y))/n
lambda_min <- lambda_max/(2e3)
k <- (0:(num_lambda - 1))/num_lambda
lambda <- lambda_max * (lambda_min/lambda_max)^k
mod <- glmnet::glmnet(X_aug, y, family = "gaussian", lambda= lambda, intercept=T,
                      standardize=F, standardize.response=F)
first_entry <- function(vect) {
  index <- ifelse(sum(abs(vect) > 0) == 0, 0, min(which(abs(vect) > 0)))
  return(index)
}
lam_indices <- apply(mod$beta, 1, function(x) first_entry(x))
lambda <- c(0, mod$lambda)
lam_vals <- lambda[lam_indices+1]*n
plot_dat <- data.frame(X = lam_vals[1:100], Xk = lam_vals[101:200], 
                       true = c(rep("Non-null features", 20), rep("Null features",80)))


ggplot(plot_dat, aes(x = X, y= Xk, color = true, shape = true)) + 
  geom_point() + scale_color_manual(values = c("red","black")) + 
  scale_shape_manual(values = c("square","circle")) +
   theme(legend.position = c(0.2, 0.9), legend.title = element_blank()) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab(TeX('Value of $\\lambda$ when $X_j$ enters the model')) + 
  ylab(TeX('Value of $\\lambda$ when $\\widetilde{X}_j$ enters the model')) +
  geom_segment(aes(x = 0, y = 1.5, xend = 1.5, yend = 1.5)) +
  geom_segment(aes(x = 1.5, y = 0, xend = 1.5, yend = 1.5)) +
  annotate(geom="text", x=2.5, y=1.8, label="Denominator") + 
  annotate(geom="text", x=2.5, y=1.6, label=TeX("($W_j \\geq t$)")) +
  annotate(geom="text", x=.5, y=2.2, label="Numerator") + 
  annotate(geom="text", x=.5, y=2.0, label=TeX("($W_j \\leq -t$)")) + 
  ggtitle('Estimated FDP at threshold t = 1.5') + 
  scale_x_continuous(limits = c(0,3),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,3),expand = c(0,0))
ggsave("simulations/figs/fig1.png")
