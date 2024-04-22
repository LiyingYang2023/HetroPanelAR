###### Pesaran and Yang (2024) "Hetro. AR panels"
rm(list=ls())

## Please change the path to the directory where the simulation outcomes are stores. 
# setwd('~/Downloads') 

## Please install the following packages. 
# list.of.packages <- c("parallel","parallelly")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(parallel) # mclapply
library(parallelly) # detect cores


#### Compare FDAC and MSW estimators
##############################################################################
#### Parameter values
##############################################################################
rep = 1000 # number of replications in each loop
reps = 1:rep
Nlist = c(100,1000); 
Tlist = c(4,6,10); 
Nmax = 5000; Tmax = 10
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])

# parameters of AR(1) coefficients
NumPara = 3
p3 = c(0.5, 0.5, 0.5, 0.5, 0) 
p4 = c(0, 0, 0.8, 1, 0)
p5 = c(0, 0, 0.85, 0.95, 0)
mphi_l = c(0.4, 0.5, p3[3]*p5[3]+p4[3]*(1-p5[3]), p3[4]*p5[4]+p4[4]*(1-p5[4]), 0.5) # mean
vphi_l = c(p3[1]^2/3, p3[2]^2/3, (p3[3]^2*p5[3]+p4[3]^2*(1-p5[3]))-mphi_l[3]^2, (p3[4]^2*p5[4]+p4[4]^2*(1-p5[4]))-mphi_l[4]^2, 0) # variance
pm = rbind(mphi_l,vphi_l,p3,p4,p5)

# parameters of initial processes
init_l= 1

# error processes
ng = matrix(c(1,0,1,0,0,0,1,1),4,2)
mcvs = cbind(0.6, 0.2) # parameters for GARCH effects
##############################################################################


##############################################################################
#### Data generating process
##############################################################################
# Generate M_i - the period of initial values
GenM = function(N,T,case,init,m0) {
  if (init == 1) {
    mvm = matrix(m0, N, 1)
    if (case == 4) {
      mv = c(m0,1)
      pmv = c(0.95, 0.05) # M_i = 2 for units with phi_i = 1
      mvm = rep(mv[1],N*pmv[1])
      for (j in 2:length(mv)) {
        mvm = c(mvm,rep(mv[j],N*pmv[j]))
      }
      mvm = matrix(t(matrix(mvm,N/100,100)),N,1)
    }
  }
  if (init == 2) {
    mv = c(10,15,20,25,30,2)
    pmv = c(0.2,0.25,0.2,0.1,0.1,0.15)
    mvm = rep(mv[1],N*pmv[1])
    for (j in 2:length(mv)) {
      mvm = c(mvm,rep(mv[j],N*pmv[j]))
    }
    mvm = matrix(t(matrix(mvm,N/100,100)),N,1)
  }
  
  Mn = matrix(1,N,T)
  for (i in 1:N) {
    Mn[i,1:(m0+1-mvm[i])] = rep(0,m0+1-mvm[i]) 
  }
  Mn
}

# Generate random coefficients
GenRandomPara_exp = function(N,NumPara,case,ve,b,kappa) {
  para = array(0,dim=c(N,NumPara))
  ## AR(1) coefficients
  if (case==1) { 
    mphi = 0.4; a = 0.5 # 0.4 + uniform (-0.5, 0.5)
    para[,3] = mphi + runif(N, min=-a, max=a) 
  }
  if (case==2) {
    mphi = 0.5; a = 0.5 # 0.5 + uniform (-0.5, 0.5)
    para[,3] = mphi + runif(N, min=-a, max=a) 
  }
  if (case==3) {
    lb = 0.5; ub = 0.8; pi = 0.85 # categorical
    para[,3] = matrix(t(matrix(c(rep(lb,pi*N),rep(ub,(1-pi)*N)),N/100,100)),N,1)
  }
  if (case==4) {
    lb = 0.5; ub = 1; pi = 0.95 # categorical
    para[,3] = matrix(t(matrix(c(rep(lb,pi*N),rep(ub,(1-pi)*N)),N/100,100)),N,1)
  }
  if (case==5) {
    mphi = 0.5 # homogeneous phi_i = phi_0 = 0.5 for all i
    para[,3] = rep(mphi,N)
  }
  ## individual-specific means: mu_i
  para[,2] = para[,3] + rnorm(N,0,1) 
  ## initial values: y_{i, -M_{i}}
  para[,1] = para[,2] + rnorm(N, b, kappa*ve)  
  para
}

# Generate the error processes
GenError_exp = function(N,T,ve,GAUSSIAN,GARCH,Mn){
  u = array(0,dim=c(N,T))
  e = array(0,dim=c(N,T))
  h = array(0,dim=c(N,T))
  if (GAUSSIAN == TRUE ){
    e[,1:T] = matrix(rnorm(N*T),N,T)
  } else {
    e[,1:T] = matrix((rchisq(N*T,2,ncp=0)-2)/2,N,T)
  }
  if (GARCH == FALSE) {
    h[,1:T] = matrix(sqrt(ve),N,T)
    u[,2:T] = Mn[,2:T]*h[,2:T]*e[,2:T] 
  } else {
    cvs = t(matrix(mcvs,2,N)) # coefficients of GARCH effects
    h[,1] = sqrt(ve) # stationary initial values for GARCH(1,1)
    for(t in 2:T){
      h[,t] = (1-Mn[,t])*sqrt(ve) + Mn[,t]*sqrt(ve-ve*rowSums(cvs) + cvs[,1]*h[,t-1]^2 + cvs[,2]*((h[,t-1]*e[,t-1])^2))
      u[,t] = Mn[,t]*h[,t]*e[,t]
    }
  }
  return(u[,2:T])
}

## DGP with alpha_i = mu_i*(1-phi_i)
GenRandomY_exp = function(N,Tmax,para,case,init,m0,ve,GAUSSIAN,GARCH){
  T = Tmax + m0 + 1
  Y = array(0,dim=c(N,T))
  Y[1:N,1] = para[1:N,1]
  Mn = GenM(N,T,case,init,m0)
  errors = GenError_exp(N,T,ve,GAUSSIAN,GARCH,Mn) # N*(T-1) matrix
  for(t in 1:(T-1)){
    ### y_{it} = mu_i*(1-phi_i) + beta_i*y_{i,t} + u_{it} with y_{i,-M_i} = 0 for all phi_i<1 and y_{i,-2} = mu_i for phi_i=1 
    Y[,(t+1)] = (1-Mn[,(t+1)])*para[,1] + Mn[,(t+1)]*(para[,2]*(1-para[,3]) + para[,3] * Y[,t] + errors[,t])
  }
  Y[,(m0+1+1):(m0+1+Tmax)]
}
##############################################################################


##############################################################################
####### MSW estimator by Mavroeidis, Sasaki and Welch (2015)
##############################################################################
### Function Kernel
K = function(u){
  # 0.75 * (1-u^2) * (u > -1) * (u < 1)
  dnorm(u) * (u != 0)
}

### Function Square Integral of Kernel
K2 = function(res){
  list = (-1*3*res):(3*res) / res
  sum( K(list)^2 ) / (res)
}

### Function Row Product
RowProd = function(U){
  nc = ncol(U)
  ones = as.matrix(array(1,nc))
  u = log(U) %*% ones
  exp(u)
}

### Function Multivariate Normal Random Number Generator
rmultnorm = function(n, mu, vmat, tol=1e-07){
  p = ncol(vmat)
  if (length(mu)!=p)
    stop("mu vector is the wrong length")
  if (max(abs(vmat - t(vmat))) > tol)
    stop("vmat not symmetric")
  vs = svd(vmat)
  vsqrt = t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
  ans = matrix(rnorm(n * p), nrow=n) %*% vsqrt
  ans = sweep(ans, 2, mu, "+")
  dimnames(ans) = list(NULL, dimnames(vmat)[[2]])
  return(ans)
}

### Function Multivariate Normal Density
dmultnorm = function(X,M,V){
  dim = ncol(X)
  d = exp( -0.5 * ( t(X-M)%*%solve(V)%*%(X-M) ) )
  d / ( (2*3.1415)^dim * abs(det(V)) )^.5
}


### Function Estimate Parameters for Time Series (Takes ROW Vector of Series Y)
GetParaHat = function(Y){
  Y = t( as.matrix(Y) )
  T = ncol(as.matrix(Y))
  Xmat = t( rbind(1, Y[1:(T-1)]) )
  Ymat = t( Y[2:T] )
  
  alpha_beta = solve(t(Xmat)%*%Xmat) %*% t(Xmat)%*%t(Ymat)
  sigma = (t(t(Ymat)-Xmat%*%alpha_beta) %*% (t(Ymat)-Xmat%*%alpha_beta) / T )^.5
  rbind(alpha_beta,log(sigma))
}

### Function Generate A Sample of (alpha hat, beta hat, log sigma hat)
GenSample = function(N,T,NumPara,M,S,r){
  para = GenRandomPara(N,NumPara,M,S,r)
  Y = GenRandomY(N,T,para)
  
  para_hat = array(0,dim=c(N,NumPara-1))
  for(i in 1:N){
    para_hat[i,1:(NumPara-1)] = t( GetParaHat(Y[i,1:T]) )
  }
  cbind(para[1:N,1] , para_hat)
}

### Function List of Parameters over Which to Integrate
GetList = function(NumPara, MinInt, MaxInt, ResInt){
  length = ResInt^NumPara
  list = array(0,dim=c(length,NumPara))
  idx_list = array(0,dim=c(length,NumPara))
  count = 1:(length) - 1
  for(i in 1:NumPara){
    idx_list[1:length,i] = count-floor(count/ResInt)*ResInt
    count = floor(count/ResInt)
  }
  idx_list = idx_list + 1
  # single_list = (1:ResInt - 0.5) * (MaxInt-MinInt)/ResInt + MinInt
  for(i in 1:NumPara){
    single_list = (1:ResInt - 0.5) * (MaxInt[i]-MinInt[i])/ResInt + MinInt[i]
    list[1:length,i] = single_list[idx_list[1:length,i]]
  }
  list
}

### Function Calculate Integral Operator (para_list - N matrix)
IntOperator = function(N,T,int_para_list,paneldata,sigma){
  list_length = nrow(int_para_list)
  alpha_list = int_para_list[1:list_length,1]
  beta_list = int_para_list[1:list_length,2]
  OperF = array(1,dim=c(list_length,N))
  
  for(i in 1:N){
    for(j in 2:T){
      OperF[1:list_length,i] = OperF[1:list_length,i] * dnorm( (paneldata[i,j] - alpha_list - beta_list * paneldata[i,j-1]) / sigma ) / sigma
    }
  }
  
  OperF
}

### Function Calculate Likelihood Integrand (para_list vector)
Integrand = function(int_para_list,MSR_hat,NumPara){
  # NumPara = 3
  list_length = nrow(int_para_list)
  M_hat = MSR_hat[1:NumPara]
  S_hat = MSR_hat[(NumPara+1):(2*NumPara)]
  R_hat = array(1,dim=c(NumPara,NumPara))
  V_hat = S_hat %*% t(S_hat)
  idx = 1
  for(i in 2:NumPara){
    for(j in 1:(i-1)){
      R_hat[i,j] = MSR_hat[2*NumPara+idx]
      R_hat[j,i] = R_hat[i,j]
      idx = idx + 1
    }
  }
  V_hat = V_hat * R_hat
  
  d_list = array(0,list_length)
  for(i in 1:list_length){
    d_list[i] = dmultnorm(as.matrix(int_para_list[i,1:NumPara]),M_hat,V_hat)
  }
  d_list
}

### Function Calculate Negative Log Likelihood
NegativeLogLikelihood = function(MSR_hat,N,T,int_para_list,paneldata,Y1,y1,h){
  NumPara = ncol(int_para_list)
  Likelihood_Integrand = Integrand(int_para_list,MSR_hat,NumPara)
  ### Below is A Non_Weighted Likelihoods across i ###
  OperF = IntOperator(N,T,int_para_list,paneldata,MSR_hat[6])
  OperF = t(OperF)
  Likelihood = OperF %*% as.matrix(Likelihood_Integrand)
  ### Below is to Avoid Infinite Log Likelihood
  # Likelihood = Likelihood + 0.0000000001
  ### Below is A Weighted Sum of Log Likelihoods across i ###
  pentalty = 100000000*sum((MSR_hat[3:4]<0)*(MSR_hat[3:4]^2)) + 0.000001*sum(MSR_hat[3:4]^(-2)) + 0.001*MSR_hat[6]^(-2) + 100000000*(MSR_hat[5]^2>1)*((MSR_hat[5]-1)^2) + .01*MSR_hat[5]^4## PENALTY ## IN THIS NEW CODE 0.0005*sum(MSR_hat[3:4]^(-2))
 -1 * t(log(Likelihood)) %*% as.matrix( K((Y1- y1)/h) ) / N / h + pentalty
}


### Function Variance of Derivative of Log Likelihood
VarDLogLikelihood = function(MSR_hat,N,T,int_para_list,paneldata,Y1,y1,h){
  NumPara = ncol(int_para_list)
  Likelihood_Integrand = Integrand(int_para_list,MSR_hat,NumPara)
  ### Below is A Non_Weighted Likelihoods across i ###
  OperF = IntOperator(N,T,int_para_list,paneldata,MSR_hat[6])
  OperF = t(OperF)
  Likelihood = OperF %*% as.matrix(Likelihood_Integrand)
  
  ddLogLikelihood = array(0,dim=c(N,6))
  
  delta = 0.0001
  for(i in 1:6){
    e = array(0,6)
    e[i] = delta
    DLikelihood_Integrand = Integrand(int_para_list,as.vector(MSR_hat+e),NumPara)
    OperF = IntOperator(N,T,int_para_list,paneldata,as.vector(MSR_hat+e)[6])
    OperF = t(OperF)
    DLikelihood = OperF %*% as.matrix(DLikelihood_Integrand)
    ddLogLikelihood[,i] = ( log(DLikelihood) - log(Likelihood) ) / delta
  }
  
  ### Conditional Mean of Derivative of Log Likelihood
  CMean = t(ddLogLikelihood) %*% as.matrix( K((Y1- y1)/h) ) / N / h
  CMean = t( array(CMean,dim=c(6,N)) )
  
  ### Conditional Mean of Derivative of Log Likelihood
  t(ddLogLikelihood - CMean) %*% diag( as.vector( K((Y1- y1)/h) ) ) %*% (ddLogLikelihood - CMean) / N / h
}

#### Parameter values for the MSW estimator
ResInt = 10
NumSim = 1000  # Number of simulations for integral operator

### Multivariate Distribution of TRUE Random Parameters 
my1 = 1     #Y_1
ma  = 0.5   #Alpha
mb  = 0.5   #Beta
M = c(my1,ma,mb)

sy1 = 0.5   #Y_1
sa  = 0.5   #Alpha
sb  = 0.5   #Beta
S = c(sy1,sa,sb)

corr = 0.5

sigma = 0.5

### Initial Values of Estimates (conditional on Y_1 = y1) 
M_hat = c(0.5,0.5) # NumPara=3 Conditional on Y_1=y1
S_hat = c(0.5,0.5) # NumPara=3 Conditional on Y_1=y1
R_hat = c(0.5)     # NumPara*(NumPara)/2 = 1
Sigma_hat = c(0.5)
MSR_hat = rbind(as.matrix(M_hat),as.matrix(S_hat),as.matrix(R_hat),as.matrix(Sigma_hat))

Kernel_mv = function(paneldata){
  Nk = max(dim(paneldata)) #  Large N fixed short T panel data
  Tk = min(dim(paneldata))
  h = 1.06 * sy1 * Nk^(-1/5)
  NumPara = NumPara #3
  MinInt = c(-1,-1)    # para_means - 2* para_stdvs
  MaxInt = c( 2, 2)    # para_means + 2* para_stdvs
  ResInt = 10
  least_sum_bias2 = 999999
  
  Y1_all = paneldata[,1] 
  idx = abs( Y1_all - mean(Y1_all) ) < 5 * var(Y1_all)^.5
  Y1 = matrix(Y1_all[idx],length(Y1_all[idx]),1)
  
  #########
  ###### Local Kernel estimators
  ########
  KernelLocal = function(y1){
    int_para_list = GetList(NumPara-1, MinInt, MaxInt, ResInt)
    results <- nlm(NegativeLogLikelihood,c(.5+.5*(y1-1),.5+.5*(y1-1),.5,.5,.5,.5),Nk,Tk,int_para_list,paneldata,Y1_all,y1,h,iterlim=5)
    results$estimate 
  }
  
  temporary_estimates = t(sapply(Y1,KernelLocal)) # Rep*6 matrix of estimators conditional on y1 values
  
  trim_indicator = temporary_estimates[,1] < 99999
  
  mc_ests = apply(temporary_estimates[trim_indicator,],2,mean) # mean of each column
  mc_ests_sq = apply(temporary_estimates[trim_indicator,]^2,2,mean)
  mc_stdev_ests = apply(temporary_estimates[trim_indicator,],2,sd)/sqrt(Nk) # non-parametric estimators of stds
  mc_stdev_ests_sq = apply(temporary_estimates[trim_indicator,]^2,2,sd)/sqrt(Nk)
  
  mc_ests_pop_vars= cbind(mc_ests_sq[1] - mc_ests[1]^2 + mc_ests_sq[3], mc_ests_sq[2] - mc_ests[2]^2 + mc_ests_sq[4])
  Valpha = cbind( cov(temporary_estimates[trim_indicator,1]^2,temporary_estimates[trim_indicator,1]^2) , cov(temporary_estimates[trim_indicator,1]^2,temporary_estimates[trim_indicator,1]) , cov(temporary_estimates[trim_indicator,1]^2,temporary_estimates[trim_indicator,3]) )
  Valpha = rbind(Valpha, cbind( cov(temporary_estimates[trim_indicator,1],temporary_estimates[trim_indicator,1]^2) , cov(temporary_estimates[trim_indicator,1],temporary_estimates[trim_indicator,1]) , cov(temporary_estimates[trim_indicator,1],temporary_estimates[trim_indicator,3]) ) )
  Valpha = rbind(Valpha, cbind( cov(temporary_estimates[trim_indicator,3],temporary_estimates[trim_indicator,1]^2) , cov(temporary_estimates[trim_indicator,3],temporary_estimates[trim_indicator,1]) , cov(temporary_estimates[trim_indicator,3],temporary_estimates[trim_indicator,3]) ) )
  Vbeta = cbind( cov(temporary_estimates[trim_indicator,2]^2,temporary_estimates[trim_indicator,2]^2) , cov(temporary_estimates[trim_indicator,2]^2,temporary_estimates[trim_indicator,2]) , cov(temporary_estimates[trim_indicator,2]^2,temporary_estimates[trim_indicator,4]) )
  Vbeta = rbind(Vbeta , cbind( cov(temporary_estimates[trim_indicator,2],temporary_estimates[trim_indicator,2]^2) , cov(temporary_estimates[trim_indicator,2],temporary_estimates[trim_indicator,2]) , cov(temporary_estimates[trim_indicator,2],temporary_estimates[trim_indicator,4]) ) )
  Vbeta = rbind(Vbeta , cbind( cov(temporary_estimates[trim_indicator,4],temporary_estimates[trim_indicator,2]^2) , cov(temporary_estimates[trim_indicator,4],temporary_estimates[trim_indicator,2]) , cov(temporary_estimates[trim_indicator,4],temporary_estimates[trim_indicator,4]) ) )
  mc_stdev_pop_vars =cbind( cbind(1, -2*mc_ests[1], 1) %*% Valpha %*% rbind(1, -2*mc_ests[1], 1), cbind(1, -2*mc_ests[2], 1) %*% Vbeta %*% rbind(1, -2*mc_ests[2], 1))
  mc_stdev_pop_vars= ( mc_stdev_pop_vars/Nk )^.5
  
  # report the estimates and the corresponding estimated standard errors
  matrix(cbind(mc_ests[2],mc_stdev_ests[2], mc_ests_pop_vars[2],mc_stdev_pop_vars[2]),nrow = 1,ncol = 4)
}
################################################################################

##############################################################################
####### FDAC and HetroGMM estimators of moments of by Pesaran and Yang (2023)
##############################################################################
inv = function(m) class(try(solve(m),silent=T))[1]=="matrix"

Moment_mv = function(paneldata){
  T = dim(paneldata)[2]; N = dim(paneldata)[1]
  y = matrix(t(paneldata),nrow=T,ncol=N)
  Tm = min(dim(y))
  tm = min(dim(y)) -1
  n = max(dim(y))
  dy = y[2:Tm,] - y[1:(Tm-1),]  # first difference of panel data
  
  ######## MM: Individual Moments average over T : (T-h-1)*N matrix 
  lagdyy = 1*(tm>2)+1*(tm>3)+1*(tm>4)+1 # order of h for \Delta y_{it}\Delta y_{i,t-h}
  dyyi1 = matrix(, nrow = 5, ncol = n) # moments of \Delta y_{it}\Delta y_{i,t-h}, h=0,1,...,lagdyy
  dyyi2 = matrix(, nrow = 3, ncol = n) # T-3 h=0,1,2 
  dyyi3 = matrix(, nrow = 4, ncol = n) # T-4 h=0,1,2,3  
  dyyi4 = matrix(, nrow = 5, ncol = n) # T-5 h=0,1,2,3,4 
  
  for (ll in 0:lagdyy){
    dyyi1[(ll+1),] = matrix(t(dy[(1+ll):tm,]*dy[1:(tm-ll),]),nrow=n,ncol=(tm-ll)) %*% matrix(1,nrow = tm-ll,ncol=1) / (tm-ll)# average over t-h-1 
  }
  
  if (tm>=3) {
    dyyi2 = rbind(matrix(matrix(t(dy[3:tm,]^2),nrow=n,ncol=(tm-2)) %*% matrix(1,nrow=tm-2,ncol=1)/(tm-2),nrow=1,ncol=n),
                  matrix(matrix(t(dy[3:tm,]*dy[2:(tm-1),]),nrow=n,ncol=(tm-2)) %*% matrix(1,nrow=tm-2,ncol=1)/(tm-2),nrow=1,ncol=n),
                  matrix(matrix(t(dy[3:tm,]*dy[1:(tm-2),]),nrow=n,ncol=(tm-2)) %*% matrix(1,nrow=tm-2,ncol=1)/(tm-2),nrow=1,ncol=n)) # MM b
    hi1 = matrix(dy[3:tm,]^(2) + dy[3:tm,]*dy[2:(tm-1),],nrow=(tm-2),ncol=n) # GMM (T-3)*n vector
    gi1 = matrix(dy[3:tm,]^(2) + 2*dy[3:tm,]*dy[2:(tm-1),] + dy[3:tm,]*dy[1:(tm-2),],nrow=(tm-2),ncol=n) 
    
    ####### MM Average across i
    dyy1 = dyyi1 %*% matrix(1,nrow=n,ncol=1)/n # h=0,1,2,...,lagdyy
    dyy2 = dyyi2 %*% matrix(1,nrow=n,ncol=1)/n # h=0,1,2
    ####### MM Version a: dyy1 with different T-h-1
    mm1 = (dyy1[1]+2*dyy1[2]+dyy1[3])/(dyy1[1]+dyy1[2]) # Equation (15) 
    mms1 = sum((dyyi1[1,]+2*dyyi1[2,]+dyyi1[3,]-mm1*(dyyi1[1,]+dyyi1[2,]))^2)/n # average across i
    std_mm1 = sqrt((dyy1[1]+dyy1[2])^(-2)*mms1/n) # ((mean(dyyi1[1,]+dyyi1[2,]))^(2)*(mms1^(-1)))^(-1)/n   
    
    
    ####### GMM 
    hnt1 = hi1 %*% matrix(1,nrow=n,ncol=1) / n  # (T-3)*1 average over n rowMeans(hi1)
    gnt1 = gi1 %*% matrix(1,nrow=n,ncol=1) / n  # average over n rowMeans(gi1) 
    gmms1 = (gi1-mm1*hi1) %*% t(gi1-mm1*hi1) /n  # the optimal weight (T-3)*(T-3) matrix for the mean Equation (20)
    if (inv(gmms1) == 1) {
      gmms1 = gmms1
    } else {
      gmms1 = diag(tm-2)/n
    }
    gmm1 =  solve(t(hnt1) %*% solve(gmms1) %*% hnt1) * (t(hnt1) %*% solve(gmms1) %*% gnt1) # 1st moment: Equation (22)
    std_gmm1 = sqrt(solve(t(hnt1) %*% solve(gmms1) %*% hnt1) /n) #Asymptotic variance
  }
  
  if (tm>=4) {
    svmi = matrix(, nrow = 2^2 , ncol=n) # sample variance of moment conditions
    svgi = matrix(, nrow = (tm-2+tm-3)^2 , ncol=n) # sample variance of joint moment conditions: 1st, 2nd
    
    dyyi3 = rbind(matrix(matrix(t(dy[4:tm,]^2),nrow=n,ncol=(tm-3)) %*% matrix(1,nrow=tm-3,ncol=1)/(tm-3),nrow=1,ncol=n), 
                  matrix(matrix(t(dy[4:tm,]*dy[3:(tm-1),]),nrow=n,ncol=(tm-3)) %*% matrix(1,nrow=tm-3,ncol=1)/(tm-3),nrow=1,ncol=n),
                  matrix(matrix(t(dy[4:tm,]*dy[2:(tm-2),]),nrow=n,ncol=(tm-3)) %*% matrix(1,nrow=tm-3,ncol=1)/(tm-3),nrow=1,ncol=n),
                  matrix(matrix(t(dy[4:tm,]*dy[1:(tm-3),]),nrow=n,ncol=(tm-3)) %*% matrix(1,nrow=tm-3,ncol=1)/(tm-3),nrow=1,ncol=n))
    hi2 = matrix(dy[4:tm,]^(2) + dy[4:tm,]*dy[3:(tm-1),],nrow=(tm-3),ncol=n) #  (T-4)*n vector
    gi2 = matrix(dy[4:tm,]^(2) + 2*dy[4:tm,]*dy[3:(tm-1),] + 2*dy[4:tm,]*dy[2:(tm-2),] + dy[4:tm,]*dy[1:(tm-3),],nrow=(tm-3),ncol=n) 
    
    dyy3 = dyyi3 %*% matrix(1,nrow=n,ncol=1)/ n # h=0,1,2,3
    ### MM: 2nd
    mm2 = (dyy1[1]+2*dyy1[2]+2*dyy1[3]+dyy1[4])/(dyy1[1]+dyy1[2]) # Equation (16) 
    mms2 = mean((dyyi1[1,]+2*dyyi1[2,]+2*dyyi1[3,]+dyyi1[4,]-mm2*(dyyi1[1,]+dyyi1[2,]))^2)
    std_mm2 = sqrt((dyy1[1]+dyy1[2])^(-2)*mms2 /n)  #1-(rho[1]+2*mm1*rho[1]+mm2*rho[1])/(1-rho[1]) - mm2 
    mmv = mm2-mm1^2 # plug-in estimator of variance
    
    mc = rbind(mm1 * (dyyi1[1,]+dyyi1[2,]) - (dyyi1[1,]+2*dyyi1[2,]+dyyi1[3,]), mm2*(dyyi1[1,]+dyyi1[2,]) - (dyyi1[1,]+2*dyyi1[2,]+2*dyyi1[3,]+dyyi1[4,])) # n*2 matrix
    index = 1
    for (i in 1:2){
      for (j in 1:2) {
        svmi[index,] = mc[i,] * mc[j,]
        index = index +1
      }
    }
    rm(mc)
    jmms2 = matrix(c(rowMeans(svmi)),nrow = 2 ,ncol = 2) # with the mmb plugged in  
    dmm2 = matrix(c(-2*mm1,1),nrow = 2,ncol = 1)
    jhm2 = matrix(c((dyy1[1]+dyy1[2]),0,0,(dyy1[1]+dyy1[2])),nrow=2,ncol=2)
    
    if (inv(jhm2) == 1) {
      jhm2 = jhm2} else {
        jhm2 = diag(2)
      }
    std_mmv = sqrt(t(dmm2) %*% solve(jhm2) %*% (jmms2) %*% solve(jhm2) %*% dmm2 /n  )
    
    ### GMM : 2nd
    hnt2 = rowMeans(hi2)
    gnt2 = rowMeans(gi2)
    gmms2 = (gi2-mm2*hi2) %*% t(gi2-mm2*hi2) /n  # the optimal weight (T-4)*(T-4) matrix for the mean Equation (31)
    if (inv(gmms2) == 1) {
      gmms2 = gmms2} else {
        gmms2 = diag(tm-3)
      }
    gmm2 =  solve(t(hnt2) %*% solve(gmms2) %*% hnt2) * (t(hnt2) %*% solve(gmms2) %*% gnt2) # 2nd moment: Equation (29)
    gmmv = gmm2 - gmm1^2
    zeroh1 = matrix(rep(0,length(hi1[,i])), nrow=length(hi1[,i]),ncol=1)
    zeroh2 = matrix(rep(0,length(hi2[,i])), nrow=length(hi2[,i]),ncol=1)
    zerom1 = matrix(rep(0, (tm-2)*(tm-3)),nrow=tm-2,ncol=tm-3)
    zerom2 = matrix(rep(0, (tm-2)*(tm-3)),nrow=tm-3,ncol=tm-2)
    
    mc = rbind(kronecker(gmm1,hi1)-gi1, kronecker(gmm2,hi2)-gi2) # (T-3+T-4)*2 matrix
    index = 1
    for (i in 1:(2*tm-5)){
      for (j in 1:(2*tm-5)) {
        svgi[index,] = mc[i,] * mc[j,]
        index = index +1
      }
    }
    rm(mc)
    jgmms2 = matrix(c(rowMeans(svgi)),nrow = 2*tm-5 ,ncol = 2*tm-5) # covariance matrix of joint moments
    jhg2 =  rbind(cbind(hnt1,zeroh1),cbind(zeroh2,hnt2)) # （2*tm-5) times 2 matrix
    
    gmma2 = rbind(cbind(solve(gmms1),zerom1),cbind(zerom2,solve(gmms2))) # weight matrix
    gmmcv2 = solve(t(jhg2) %*% gmma2 %*% jhg2) %*% t(jhg2) %*% gmma2 %*% jgmms2 %*% t(gmma2) %*% jhg2 %*% solve(t(jhg2) %*% t(gmma2) %*% jhg2)
    dgmm2 = matrix(c(-2*gmm1,1),nrow = 2,ncol = 1)
    
    std_gmm2 = sqrt((t(hnt2) %*% solve(gmms2) %*% hnt2)^(-1)/n)
    std_gmmv =  sqrt(t(dgmm2) %*% gmmcv2  %*% dgmm2 /n)
    
  }
  
  if (tm>=5) {
    dyyi4 = rbind(matrix(matrix(t(dy[5:tm,]^2),nrow=n,ncol=(tm-4)) %*% matrix(1,nrow=tm-4,ncol=1)/(tm-4),nrow=1,ncol=n),
                  matrix(matrix(t(dy[5:tm,]*dy[4:(tm-1),]),nrow=n,ncol=(tm-4)) %*% matrix(1,nrow=tm-4,ncol=1)/(tm-4),nrow=1,ncol=n),
                  matrix(matrix(t(dy[5:tm,]*dy[3:(tm-2),]),nrow=n,ncol=(tm-4)) %*% matrix(1,nrow=tm-4,ncol=1)/(tm-4),nrow=1,ncol=n),
                  matrix(matrix(t(dy[5:tm,]*dy[2:(tm-3),]),nrow=n,ncol=(tm-4)) %*% matrix(1,nrow=tm-4,ncol=1)/(tm-4),nrow=1,ncol=n),
                  matrix(matrix(t(dy[5:tm,]*dy[1:(tm-4),]),nrow=n,ncol=(tm-4)) %*% matrix(1,nrow=tm-4,ncol=1)/(tm-4),nrow=1,ncol=n))
    hi3 = matrix(dy[5:tm,]^(2) + dy[5:tm,]*dy[4:(tm-1),],nrow=(tm-4),ncol=n) #  (T-5)*n vector
    gi3 = matrix(dy[5:tm,]^(2) + 2*dy[5:tm,]*dy[4:(tm-1),] + 2*dy[5:tm,]*dy[3:(tm-2),] + 2*dy[5:tm,]*dy[2:(tm-3),] + dy[5:tm,]*dy[1:(tm-4),],nrow=(tm-4),ncol=n) # 
    
    dyy4 = rowMeans(dyyi4) # h=0,1,2,3,4
    ### MM
    mm3 = (dyy1[1]+2*dyy1[2]+2*dyy1[3]+2*dyy1[4]+dyy1[5])/(dyy1[1]+dyy1[2]) # 
    mms3 = mean((dyyi1[1,]+2*dyyi1[2,]+2*dyyi1[3,]+2*dyyi1[4,]+dyyi1[5,]-mm3*(dyyi1[1,]+dyyi1[2,]))^2)
    std_mm3 = sqrt((dyy1[1]+dyy1[2])^(-2)*mms3 /n)  #1-(rho[1]+2*mm1*rho[1]+mm2*rho[1])/(1-rho[1]) - mm2 
    
    
    ### GMM: 3rd
    hnt3 = rowMeans(hi3)
    gnt3 = rowMeans(gi3)
    gmms3 = (gi3-mm3*hi3) %*% t(gi3-mm3*hi3) /n  # the optimal weight (T-5)*(T-5) matrix for the mean Equation (31)
    if (inv(gmms3) == 1) {
      gmms3 = gmms3} else {
        gmms3 = diag(tm-4)
      }
    gmm3 =  solve(t(hnt3) %*% solve(gmms3) %*% hnt3) * (t(hnt3) %*% solve(gmms3) %*% gnt3) # 2nd moment: Equation (29)
    std_gmm3 = sqrt((t(hnt3) %*% solve(gmms3) %*% hnt3)^(-1)/n)
  }
  
  if (tm==3){
    return(cbind(mm1,std_mm1,gmm1,std_gmm1))
  }
  if (tm==4) {
    return(cbind(mm1,std_mm1,gmm1,std_gmm1,
                 mm2,std_mm2,gmm2,std_gmm2,
                 mmv,std_mmv,gmmv,std_gmmv))
  }
  if (tm>4){
    return(cbind(mm1,std_mm1,gmm1,std_gmm1,
                 mm2,std_mm2,gmm2,std_gmm2,
                 mmv,std_mmv,gmmv,std_gmmv,
                 mm3,std_mm3,gmm3,std_gmm3))
  }
}
##############################################################################

MonteCarloSim = function(r,N,Tobs){
  NumPara = NumPara #3
  ResInt = ResInt #10
  NumSim = NumSim #1000  # Number of simulations for integral operator
  
  my1 = my1 #1     #Y_1
  ma  = ma #0.5   #Alpha
  mb  = mb #0.5   #Beta
  M = c(my1,ma,mb)
  
  sy1 = sy1 #0.5   #Y_1
  sa  = sa #0.5   #Alpha
  sb  = sb #0.5   #Beta
  S = c(sy1,sa,sb)
  
  corr = corr #0.5
  sigma = sigma #0.5
  panel0 = matrix(data_all[,r],Nmax,Tmax); 
  panel1 = panel0[1:N,1:Tobs]
  cbind(Kernel_mv(panel1),Moment_mv(panel1))
}


init = 1; b = 1; kappa = 2; GAUSSIAN = 1; GARCH = 0
for (m0 in c(100,1)) {
  for (case in c(1,2,5)) {
    set.seed(123987654) 
    #### Generate samples
    ##################################################################################################
    data_all = matrix(,Nmax*Tmax,rep); 
    for (r in 1:rep ) {
      ve = 0.5+ 0.5*rchisq(Nmax,1,ncp=0) # cross-sectional hetero variance of errors
      para = GenRandomPara_exp(Nmax,NumPara,case,ve,b,kappa)
      paneldata = GenRandomY_exp(Nmax,Tmax,para,case,init,m0,ve,GAUSSIAN,GARCH)
      data_all[,r] = matrix(paneldata,Nmax*Tmax,1);rm(paneldata)
    }
    rm(ve,para)
    ##################################################################################################
    for (Tobs in Tlist) {
      for (N in Nlist) {
        numCores = availableCores()
        mc_results = mclapply(reps,MonteCarloSim,N=N,Tobs=Tobs,mc.preschedule = TRUE, mc.cores = numCores-1)
        
        name = paste("mc","_c",case,"_m",m0,sep="")
        Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
        s = (Tid-1)*length(Nlist)+Nid;
        if (s<10) {
          name = paste(name,s,sep="_0")
        } else {
          name = paste(name,s,sep="_") 
        }
        assign(name,mc_results);rm(name,mc_results)
      }
    }
    rm(data_all)
    fn = paste("exp_fdac_msw_exp_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
    save.image(file=fn)
  } 
}








