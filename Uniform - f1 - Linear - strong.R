rm(list = ls())

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

NL_ME <- function(N, SD_ME, beta_0){
  
  out = foreach(i = 1:500) %dopar% {
    
    library(splines)
    library(mgcv)
    library(mfp)
    
    expit <- function(x){exp(x)/(1+exp(x))}
    
    seeds = read.table("seeds.txt")
    
    n <- N
    
    x.pred <- seq(100,200, 100/299)
    
    set.seed(seeds[i,1])
    
    X <- runif(n, 100, 200)
    X.star <- X + rnorm(n, 0, SD_ME)
    
    beta_1 <- log(1.8^(1/28.9))
    
    fx_1 <- beta_0 + beta_1*X
    
    # plot(x=X, y=fx_4, cex = 0.5)
    
    p_1 <- expit(fx_1)
    #mean(p_1)
    
    Y <- rbinom(n, 1, p_1)
    
    ### Regression ###
    
    reg.bs.X <- glm(Y ~ -1 + bs(X, degree=3, intercept=T, knots = c(median(X))), family = "binomial")
    reg.bs.X.star <- glm(Y ~ -1 + bs(X.star, degree=3, intercept=T, knots = c(median(X.star))), family = "binomial")
    
    fp.X <- mfp(Y ~ fp(X, select = 1, alpha = 1), family = "binomial")
    fp.X.star <- mfp(Y ~ fp(X.star, select = 1, alpha = 1), family = "binomial")
    
    ps.X <- gam(Y ~ s(X, bs="ps", k = 14), family = "binomial", method = "REML")
    ps.X.star <- gam(Y ~ s(X.star, bs="ps", k = 14), family = "binomial", method = "REML")
    
    ns.X <- glm(Y ~ ns(X, knots = c(median(X)), Boundary.knots = quantile(X, probs = c(0.05, 0.95))), family = "binomial")
    ns.X.star <- glm(Y ~ ns(X.star, knots = c(median(X.star)), Boundary.knots = quantile(X.star, probs = c(0.05, 0.95))), family = "binomial")
    
    pred_reg.bs.X <- predict(reg.bs.X, type = "link", newdata = data.frame(X = x.pred))
    pred_reg.bs.X.star <- predict(reg.bs.X.star, type = "link", newdata = data.frame(X.star = x.pred))
    
    pred_reg.FP.X <- predict(fp.X, type = "link", newdata = data.frame(X = x.pred))
    pred_reg.FP.X.star <- predict(fp.X.star, type = "link", newdata = data.frame(X.star = x.pred))
    
    pred_reg.ps.X <- predict(ps.X, type = "link", newdata = data.frame(X = x.pred))
    pred_reg.ps.X.star <- predict(ps.X.star, type = "link", newdata = data.frame(X.star = x.pred))

    pred_reg.ns.X <- predict(ns.X, type = "link", newdata = data.frame(X = x.pred))
    pred_reg.ns.X.star <- predict(ns.X.star, type = "link", newdata = data.frame(X.star = x.pred))
    
    rbind(pred_reg.bs.X, pred_reg.bs.X.star, 
          pred_reg.FP.X, pred_reg.FP.X.star,
          pred_reg.ps.X, pred_reg.ps.X.star,
          pred_reg.ns.X, pred_reg.ns.X.star)
    
  }
  
  return(out)
  
}

stopCluster(cl)
