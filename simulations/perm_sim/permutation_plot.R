library(prelim.knockoffs)
library(mvtnorm)
library(latex2exp)

normalize <- function(X) {
  X.scaled <- scale(X, center = F, scale = sqrt(colSums(X^2)))
  X_scaled <- X.scaled[, ]
  return(X_scaled)
}

generate_dat <- function() {
  n <- 300
  p <- 100
  mean <- rep(0,p)
  cov <- matrix(data = 0.3, nrow = p, ncol = p) + diag(x = 0.7, nrow = p)
  X <- rmvnorm(n, mean, cov)
  X <- normalize(scale(X, center = T, scale = F))
  #X <- scale(X, center= FALSE, scale = TRUE)
  #X <- scale(X)
  z <- rnorm(n, 0, 1) 
  x_30_sum <- apply(X, 1, function(x) sum(x[1:30]))
  y <- 3.5*x_30_sum + z
  return(list(X=X,y=y))
}

n <- 300
p <- 100
dat <- generate_dat()
X <- dat$X
y <- dat$y

perm <- sample(1:n ,n, replace = F)
X_perm <- X[perm,]

X_aug <- cbind(X, X_perm)

n <- nrow(X_aug)
p <- ncol(X_aug)
#X_aug <- scale(X_aug)[,] #standardize

#another way 
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
lam_vals <- lambda[lam_indices+1]*(n)
plot(lam_vals)

save(lam_vals, file = "simulations/perm_sim/lambda_vals.rda")

plot_dat <- data.frame(x = 1:200, y = lam_vals, 
                       true = c(rep("Original non-null features",30),
                                rep("Original null features",70),
                                rep("Permuted features",100)))
ggplot(data = plot_dat, aes(x = x, y = y, color = true, shape = true)) + 
  geom_point() + scale_color_manual(values = c("black", "red", "blue")) + 
  scale_shape_manual(values = c("square","circle","triangle")) +
  xlab(TeX("Index of column in the augmented matrix \\[X $X^{\\pi}$\\]")) +
  ylab(TeX("$\\lambda$ value when variable enters model")) + 
  theme_gray(base_size = 18) +
  theme(legend.position = c(0.7, 0.8), legend.title = element_blank()) 
ggsave("simulations/figs/fig2.png",
       width = 12, height = 5)
