gsm=function(mydata, mykin, nfold){
  
## negative loglikelihood
nloglik.REML.1d<- function(t, y, kin) {
  n=length(y)
  x.design=rep(1, times=n)
  beta=t[1]
  phi.11=exp(t[2])
  
  sigma.11=exp(t[3])
  
  v.phi=kin*phi.11
  v.sigma=diag(n)*sigma.11
  v=v.phi+v.sigma
  u=beta
  nloglik = 0.5*((n-2)*log(2*pi)-log(det(t(x.design)%*%x.design)) +log(det(v))+log(det(t(x.design)%*%solve(v)%*%x.design))++t(y-u)%*%solve(v)%*%(y-u))
  
  return(nloglik)
}


## henderson equation
henderson.FUN <- function(t, y, kin) {
  n <- length(y)
  x.design <- rep(1, times=n)
  
  phi.11 <- exp(t[2])
  sigma.11 <- exp(t[3])
  lambda <- phi.11/sigma.11
  
  mm <- matrix(NA, n+1, n+1)
  mm[1, 1] <- t(x.design)%*%x.design
  mm[2:(n+1), 1] <- x.design
  mm[1, 2:(n+1)] <- t(x.design)
  mm[2:(n+1), 2:(n+1)] <- diag(n) + solve(kin)/lambda
  
  est <- solve(mm)%*%c(t(x.design)%*%y, y)
  v.est <- solve(mm)*sigma.11
  
  beta_hat <- est[1]
  eta_hat <- est[-1]
  
  return(list(beta=beta_hat, eta=eta_hat, v=v.est))
}

cv.FUN=function(y0, kin0, nfold){
  n=length(y0)
  foldid <- sample(rep(1:nfold, n/nfold), replace=F)
  y.test.pred.cv=NULL
  y.test.cv=NULL
  for (ifold in 1:nfold) {
    print(ifold)
    y.test = y0[foldid==ifold]
    
    K.train=kin0[foldid!=ifold,foldid!=ifold]
    K.test=kin0[foldid==ifold,foldid==ifold]
    K.test.train=kin0[foldid==ifold,foldid!=ifold]
    
    y.train=y0[foldid!=ifold]
    t.init=c(mean(y.train), 1, 1)
    junk=optim(t.init, y=y.train, kin=K.train, fn=nloglik.REML.1d, method="Nelder-Mead", hessian=FALSE, control=list(fnscale=1, maxit=2000))
    ss.est=exp(junk$par[2])
    ee.est=exp(junk$par[3])
    v.train=ss.est*K.train+ee.est*diag(dim(K.train)[1])
    
    junk2=henderson.FUN(t=junk$par, y=y.train, kin=K.train)
    beta.hat=junk2$beta
    
    #### prediction
    v.test.train=ss.est*K.test.train
    y.test.pred=rep(1, length(y.test))*beta.hat+v.test.train%*%solve(v.train)%*%(y.train-rep(1, length(y.train))*beta.hat)
    
    y.test.cv=c(y.test.cv, y.test)
    
    y.test.pred.cv=c(y.test.pred.cv, y.test.pred)
  }
  
  xx=y.test.cv-rep(1, length(y.test.cv))*beta.hat
  xx.hat=y.test.pred.cv-rep(1, length(y.test.cv))*beta.hat
  
  predict=cov(xx, xx.hat)^2/var(xx)/var(xx.hat)
  return(list(predict=predict))
}

junk=optim(c(mean(mydata),1, 1), y=mydata, kin=mykin, fn=nloglik.REML.1d, method="Nelder-Mead", hessian=FALSE, control=list(fnscale=1, maxit=2000))
ss.est=exp(junk$par[2])
ee.est=exp(junk$par[3])
junk2=henderson.FUN(junk$par, mydata, mykin)
beta.est=junk2$beta
eta.est=junk2$eta
v.est=junk2$v
predic.est=cv.FUN(mydata, mykin, nfold)
res=list(beta.est=beta.est, ss.est=ss.est, ee.est=ee.est, eta.est=eta.est, v.est=v.est, predic.est=predic.est)
return(res)
}

  
