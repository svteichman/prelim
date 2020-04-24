set.seed(8)
p=100; n=200; k=15
X = matrix(rnorm(n*p),n)
nonzero = sample(p, k)
beta = 5.5 * (1:p %in% nonzero)
y = X %*% beta + rnorm(n)
their_res <- knockoff::knockoff.filter(X,y,knock_gen,stat.glmnet_lambdasmax,fdr=.2)
their_res2 <- knockoff::knockoff.filter(X,y,knock_gen,stat.glmnet_lambdasmax,fdr=.2)
my_res <- prelim.knockoffs::knockoff_filter(X,y,"sdp")
my_res$selected
their_res$selected
their_res2$selected