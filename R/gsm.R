gsm=function(mydata, mykin, nfold){
  
    nloglik.REML.1d<- function(t, data, kin) {
    y=c(data)
    n=length(data)
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
    
    y1=mydata$y1
    y2=mydata$y2
    foldid=mydata$id
    
for (ifold in 1:nfold) {
  y1.test = y1[foldid==ifold]; y2.test = y2[foldid==ifold]
  
  K.train=mykin[foldid!=ifold,foldid!=ifold]
  K.test=mykin[foldid==ifold,foldid==ifold]
  K.test.train=mykin[foldid==ifold,foldid!=ifold]
  
  y1.train=y1[foldid!=ifold]; y2.train=y2[foldid!=ifold]
  t1.init=c(mean(y1.train), 1, 1); t2.init=c(mean(y2.train), 1, 1)
  
  junk1=optim(t1.init, data=y1.train, kin=K.train, fn=nloglik.REML.1d, method="Nelder-Mead", hessian=FALSE, control=list(fnscale=1, maxit=2000))
  junk2=optim(t2.init, data=y2.train, kin=K.train, fn=nloglik.REML.1d, method="Nelder-Mead", hessian=FALSE, control=list(fnscale=1, maxit=2000))
  
  beta.hat=c(junk1$par[1], junk2$par[1])
  ss1.est=exp(junk1$par[2]); ss2.est=exp(junk2$par[2])
  ee1.est=exp(junk1$par[3]); ee2.est=exp(junk2$par[3])
  v1.train=ss1.est*K.train+ee1.est*diag(dim(K.train)[1])
  v2.train=ss2.est*K.train+ee2.est*diag(dim(K.train)[1])
  
  #### prediction
  v1.test.train=ss1.est*K.test.train; v2.test.train=ss2.est*K.test.train
  y1.test.pred=rep(1, length(y1.test))%*%as.matrix(beta.hat[1])+v1.test.train%*%solve(v1.train)%*%(y1.train-rep(1, length(y1.train))%*%as.matrix(beta.hat[1]))
  y2.test.pred=rep(1, length(y2.test))%*%as.matrix(beta.hat[2])+v2.test.train%*%solve(v2.train)%*%(y2.train-rep(1, length(y2.train))%*%as.matrix(beta.hat[2]))
  
  y1.test.cv=c(y1.test.cv, y1.test)
  y2.test.cv=c(y2.test.cv, y2.test)
  
  y1.test.pred.cv=c(y1.test.pred.cv, y1.test.pred)
  y2.test.pred.cv=c(y2.test.pred.cv, y2.test.pred)
}
    y.test.pred.cv=cbind(y1.test.pred.cv, y2.test.pred.cv)
    ss.cv=cov(y.test.pred.cv)
    residual=cbind((y1.test.cv-y1.test.pred.cv), (y2.test.cv-y2.test.pred.cv))
    ee.cv=cov(residual)
    h.cv=ss.cv/(ss.cv+ee.cv)
    res=list(ss.cv=ss.cv, ee.cv=ee.cv, h.cv=h.cv)
    return(res)
}


    
    
