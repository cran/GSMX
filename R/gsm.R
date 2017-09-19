gsm=function(mydata, mykin, nfold){
  
  ### Negative nloglikelihood Function
  nloglik.REML.1d <- function(t, y, kin){
    
    n <- length(y)
    x.design <- matrix(1, n, 1)
    
    phi <- t[1]
    
    sigma <- t[2]
    
    v.phi <- phi*kin
    v.sigma <- sigma*diag(n)
    v <- v.phi+v.sigma
    
    beta <- solve(t(x.design)%*%solve(v)%*%x.design)%*%(t(x.design)%*%solve(v)%*%y)
    
    nloglik = 0.5*(unlist(determinant(v))[[1]] + unlist(determinant(t(x.design)%*%solve(v)%*%x.design))[[1]] + t(y-x.design%*%beta)%*%solve(v)%*%(y-x.design%*%beta))
    
    return(nloglik)
  }
  
  
  ### Henderson Equation and Estimation
  henderson.FUN <- function(t, y, kin){
    n <- length(y)
    x.design <- matrix(1, n, 1)
    
    
    phi <- t[1]
    
    sigma <- t[2]
    
    lambda <- phi/sigma
    
    TL <- t(x.design)%*%x.design
    BL <- x.design
    TR <- t(x.design)
    BR <- diag(n) + solve(kin)/lambda
    
    v.est <- solve(cbind(rbind(TL, BL), rbind(TR, BR)))
    est <- v.est%*%matrix(c(t(x.design)%*%y, y), n+1, 1)
    
    beta_hat <- est[1, 1]
    eta_hat <- est[-1, 1]
    
    return(list(beta=beta_hat, eta=eta_hat, v=v.est))
  }
  
  
  ### Cross Validation
  cv.FUN <- function(y0, kin0, nfold){
    n <- length(y0)
    x.design <- matrix(1, n, 1)
    
    t.init <- c(var(y0)*2/3, var(y0)*1/3)
    #t.init <- c(1, 0.1)
    
    foldid <- sample(rep(1:nfold, n/nfold), replace=F)
    
    y.test.cv <- NULL
    y.test.pred.cv <- NULL
    
    for (ifold in 1:nfold){
      print(ifold)
      y.test <- y0[foldid==ifold]
      n.test <- length(y.test)
      
      x.design.test <- matrix(1, n.test, 1)
      
      
      K.train <- kin0[foldid!=ifold,foldid!=ifold]
      K.test <- kin0[foldid==ifold,foldid==ifold]
      K.test.train <- kin0[foldid==ifold,foldid!=ifold]
      
      y.train <- y0[foldid!=ifold]
      n.train <- length(y.train)
      
      x.design.train <- matrix(1, n.train, 1)
      
      junk <- optim(par=t.init, y=y.train, kin=K.train, fn=nloglik.REML.1d, method="L-BFGS-B", lower=c(0, 0), upper=c(1000, 1000), hessian=TRUE)
      
      phi.est <- junk$par[1]
      
      sigma.est <- junk$par[2]
      
      v.phi.est <- phi.est*K.train
      v.sigma.est <- sigma.est*diag(n.train)
      v.train <- v.phi.est + v.sigma.est
      
      junk2 <- henderson.FUN(t=junk$par, y=y.train, kin=K.train)
      beta.hat <- junk2$beta
      
      ### prediction ###
      y.test.pred <- x.design.test%*%beta.hat + phi.est*K.test.train%*%solve(phi.est*K.train+ sigma.est*diag(n.train))%*%(y.train-x.design.train%*%beta.hat)
      
      y.test.cv <- c(y.test.cv, y.test)
      
      y.test.pred.cv <- c(y.test.pred.cv, y.test.pred)
    }
    
    
    xx <- y.test.cv - x.design%*%beta.hat
    xx.hat <- y.test.pred.cv - x.design%*%beta.hat
    
    
    #xx <- y.test.cv
    #xx.hat <- y.test.pred.cv
    
    predict <- cov(xx, xx.hat)^2/(var(xx)*var(xx.hat))
    
    return(list(predict=predict))
  }
  
  
  n <- length(mydata)
  
  t.init <- c(var(mydata)*2/3, var(mydata)*1/3)
  #t.init <- c(1, 0.1)
  
  junk <- optim(par=t.init, y=mydata, kin=mykin, fn=nloglik.REML.1d, method="L-BFGS-B", lower=c(0, 0), upper=c(1000, 1000), hessian=TRUE)
  
  var.para <- junk$par
  
  junk2 <- henderson.FUN(junk$par, mydata, mykin)
  
  beta.est <- junk2$beta
  eta.est <- junk2$eta
  v.est <- junk2$v
  
  predic.est <- cv.FUN(mydata, mykin, nfold)
  res <- list(beta.est=beta.est, var.para=var.para, eta.est=eta.est, v.est=v.est, predic.est=predic.est)
  return(res)
  
}

  
