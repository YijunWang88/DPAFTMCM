#####  colorectal cancer clinical trial data #####
library(survival)
library(mnormt)
hos<-read.csv("hospitalization_data.csv")
hos1<-hos
for(i in 1:403){
  temp=hos[hos[,1]==i,]
  tempenum=max(temp[,2])
  hos1[i,]=temp[temp$enum==tempenum,]
}

hos1<-hos1[1:403,]
n<-dim(hos1)[1]
chemo1<-ifelse(hos1$chemo=="Treated",1,0)
chemo0<-1-chemo1
female=ifelse(hos1$sex=="Female",1,0)
male=1-female
dukes1=ifelse(hos1$dukes=="A-B",1,0)
dukes2=ifelse(hos1$dukes=="C",1,0)
dukes3=ifelse(hos1$dukes=="D",1,0)
charlson1=ifelse(hos1$charlson==0,1, 0)
charlson2=ifelse(hos1$charlson=="1月2日",1, 0)
charlson3=ifelse(hos1$charlson==3,1,0)
z<-cbind(chemo1,female,dukes2,dukes3,charlson1,charlson2)
x=as.data.frame(z)
Time<-hos1$t.stop
Status<-rep(1,n)
epsilon<-ifelse(hos1$death==0,1,2)
y=1-Status
if(class(x)=="numeric"|class(x)=="integer")xnum<-1
if(class(x)=="matrix"|class(x)=="data.frame")xnum<-dim(x)[2];x<-as.matrix(x)

#####specify hyperparameters for priors
p=xnum
a0=0;b0=100;
c0=.1;d0=.1;
g=.1;q=.1;
alpha=0.03
beta0=c(rep(0,length=xnum))
Sig_beta=diag(xnum)*1000#matrix(c(1000,0,0,1000),2,2)
gamma0=c(rep(0,length=xnum+1))
Sig_gamma=diag(xnum+1)*1000#matrix(c(1000,0,0,0,1000,0,0,0,1000),3,3)

loglike_gamma<-function(gamma,x=x,y=y){
  precision=700
  gx= as.numeric(cbind(1,x)%*%gamma)
  logq <- log(1+exp(gx))
  logq[gx>precision]<-gx[gx>precision]### to avoid inf
  f<- sum( y*gx - logq ,na.rm = TRUE)
  gx[gx>precision]<-precision
  X<-cbind(1,x)
  h <- 1/((1+exp(-gx))*(1+exp(gx)))
  H <- t(X) %*% diag( h ) %*% X
  return(loglike=f)
}

mh_gamma<- function(nsim,gamma0,Sig_gamma,prop_sigma_gamma,current.gamma){
  mh_gamma_out<-matrix(ncol=length(gamma0),nrow=nsim)
  acc.gamma<-0
  
  for(i in 1:nsim){
    cur<-loglike_gamma(current.gamma,x=x,y=y )
    prop.gamma<-rnorm(xnum+1,current.gamma,prop_sigma_gamma)
    prop<-loglike_gamma( prop.gamma,x=x,y=y )
    rate_gamma<-prop-cur+
      sum(dmnorm(prop.gamma,gamma0,Sig_gamma,log = TRUE))-
      sum(dmnorm(current.gamma,gamma0,Sig_gamma,log = TRUE))#+
    if(log(runif(1))<rate_gamma){
      current.gamma<-prop.gamma
      acc.gamma<-acc.gamma+1
    }
    mh_gamma_out[i,]<-current.gamma
    
  }
  return(list(mh_gamma_out = mh_gamma_out, acc.gamma = acc.gamma))
}

prop_sigma_gamma<-c(rep(0.001,length=xnum+1))
current.gamma<-c(rep(1,length=xnum+1))#runif(3)
nsim=5000
mh1<-mh_gamma(nsim,gamma0,Sig_gamma,
              prop_sigma_gamma,current.gamma)
mh_gamma_out<-mh1$mh_gamma_out

current.gamma<-apply(mh_gamma_out[-(1:(nsim/2)),],2,mean)
print(current.gamma)

nsim=50 #number of MCMC iterations
DPaft<-function(nsim,N=20,n=n,Time=Time,y=y,x=x){
  z=matrix(0,nrow=n,ncol=nsim)
  sim_alpha=matrix(0,ncol=nsim,nrow=1)
  sim_beta=matrix(0,ncol=nsim,nrow=xnum)
  sim_gamma=matrix(0,ncol=nsim,nrow=xnum+1)
  sim_taue=matrix(0,ncol=nsim,nrow=1)
  
  t=seq(-5,5,.01)
  denW=matrix(0,ncol=nsim,nrow=length(t))
  dent=matrix(0,ncol=nsim,nrow=length(t))
  tt=seq(0,max(Time),length=n)
  Dent=matrix(0,ncol=nsim,nrow=length(tt))
  
  #specify initial values for MCMC loops
  N=20
  muh=rep(NA,N)
  tauh=rep(NA,N)
  U=rep(NA,n)
  s=rep(NA,n)
  vh=rep(NA,N)
  ph=rep(1/N,N)
  ind=array(dim=c(N,n))
  # MCMC的第一步
  sim_alpha[1]=0.001
  sim_beta[,1]=beta0
  sim_gamma[,1]=gamma0
  muh=rnorm(N,a0,b0)
  tauh<-numeric(N)
  for(i in 1:N){
    tauh[i]=rgamma(1,c0,rate=d0)
    while(is.na(tauh[i])){
      tauh[i]=rgamma(1,c0,rate=d0)
    }
  }
  
  P=rep(1/N,N)
  ind=rmultinom(n,1,P) 
  for(i in 1:n){s[i]=seq(1,N)[ind[,i]==1]}
  nh=rowSums(ind)
  for(m in 1:(nsim-1)){
    z[,m]=log(Time)
    for (i in 1:N){
      new_tau=1/(tauh[i]*nh[i]+1/b0^2)
      new_mu=new_tau*(a0/b0^2+tauh[i]*sum((z[,m]+x%*%sim_beta[,m])[s==i]))
      muh[i]=rnorm(1,new_mu,sqrt(new_tau))
      tauh[i]=rgamma(1,c0+nh[i]/2,rate=d0+1/2*sum(((z[,m]+x%*%sim_beta[,m]-muh[i])[s==i])^2))
    }
    
    for(i in 1:n){U[i]=runif(1,max=ph[s[i]])}
    ustar=min(U)
    
    tempn=sum(nh)-cumsum(nh)
    vh[1]=rbeta(1,nh[1]+1,sim_alpha[1,m]+tempn[1])
    temp_v=1
    for (i in 2:N){
      vh[i]=rbeta(1,nh[i]+1,sim_alpha[1,m]+tempn[i])
      temp_v=temp_v*(1-vh[i-1])
      temp_alpha=max(U[s==i]/temp_v,1e-10)
      vh[i]=max(temp_alpha,vh[i]) 
      vh[i]=min(1-1e-10,vh[i])
      #if(N>1000)break
    }
    
    tempv=cumprod(c(1,1-vh))
    for(i in 1:N){ph[i]=tempv[i]*vh[i]}
    tempvv=cumsum(log(1-vh))
    sim_alpha[1,m+1]=rgamma(1,g+sum(c(1:N)),rate=q-sum(tempvv))
    
    for(i in 1:n){U[i]=runif(1,max=ph[s[i]])}
    ustar=min(U)
    NN=N
    while (sum(ph)<1-ustar){
      muh=c(muh,0)
      tauh=c(tauh,0)
      vh=c(vh,0)
      ph=c(ph,0)
      muh[N+1]=rnorm(1,a0,b0)
      tauh[N+1]=rgamma(1,c0,rate=d0)
      vh[N+1]=rbeta(1,1,sim_alpha[1,m+1])
      vh[N+1]=min(1-1e-10,vh[N+1])
      ph[N+1]=ph[N]*vh[N+1]*(1-vh[N])/vh[N]
      N=N+1
      #if(N>1000)break
    }
    
    ind=array(dim=c(N,n))
    for (i in (1:n)){ 
      ind[,i]=rmultinom(1,1,(dnorm(z[i,m],as.vector(-x[i,]%*%sim_beta[,m])+muh,1/sqrt(tauh)))*(ph>U[i]))
    }
    for (i in 1:length(t)){denW[i,m]=sum(ph*dnorm(t[i],muh,1/sqrt(tauh)))}
    
    xm=as.vector(apply(x,2,median))
    for(i in 1:length(t)){dent[i,m]=sum(ph*dnorm(t[i],as.vector(-xm%*%sim_beta[,m])+muh+muh,1/sqrt(tauh)))}
    
    for(i in 1:length(tt)){Dent[i,m]=sum(ph*pnorm(log(tt[i]),as.vector(-xm%*%sim_beta[,m])+muh,1/sqrt(tauh)),na.rm=TRUE)}
    
    for(i in 1:n){s[i]=seq(1,N)[ind[,i]==1]}
    nh=rowSums(ind)
    
    Sig=solve(solve(Sig_beta)+t(tauh[s]*x)%*%x)
    Beta=Sig%*%(solve(Sig_beta)%*%beta0+t(x)%*%(tauh[s]*(-z[,m]+muh[s])))
    
    SigE=eigen(Sig)
    SigEV=diag(SigE$values)
    SigH=SigE$vectors%*%sqrt(SigEV)
    
    sim_beta[,m+1]=Beta+SigH%*%rnorm(xnum)
    #print(m)
  }
  #mh1<-mh_gamma(nsim,gamma0,Sig_gamma,
  #              prop_sigma_gamma,current.gamma)
  #sim_gamma<-mh1$mh_gamma_out
  return(list(sim_beta=sim_beta,sim_alpha=sim_alpha,dent=dent,
              denW=denW,Dent=Dent,N=N,ind=ind,muh=muh,
              tauh=tauh,vh=vh,ph=ph))
}


nsim=5000
emlist <- DPaft(nsim=nsim,N=20,n=n,Time=Time,y=y,x=x)
beta<-apply(emlist$sim_beta[,-(1:(nsim/2))],1,mean)
mh1<-mh_gamma(nsim,gamma0,Sig_gamma,
              prop_sigma_gamma,current.gamma)
gamma<-apply(mh1$mh_gamma_out[-(1:(nsim/2)),],2,mean)
gamma<-current.gamma
para<-c(gamma,beta)
