#### FDAC and HetroGMM estimators proposed in Pesaran and Yang (2024) "Hetro. AR panels"
rm(list=ls())

## Please change the path to the directory where the simulation outcomes are stores. 
# setwd('~/Downloads') 

## Please install the following packages. 
# list.of.packages <- c("parallel","parallelly")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(parallel) # mclapply
library(parallelly) # detect cores

#### Compare FDAC and HetroGMM estimators of sigma_phi^2 = Var(phi_i)
#############################################################################
#### Parameter values
##############################################################################
NumPara = 3
Nlist = c(100,1000,2500,5000); 
Tlist = c(4,5,6,10); 
Nmax = max(Nlist); Tmax = max(Tlist)
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])

# parameters of AR(1) coefficients
p3 = c(0.5, 0.5, 0.5, 0.5, 0) 
p4 = c(0, 0, 0.8, 1, 0)
p5 = c(0, 0, 0.85, 0.95, 0)
mphi_l = c(0.4, 0.5, p3[3]*p5[3]+p4[3]*(1-p5[3]), p3[4]*p5[4]+p4[4]*(1-p5[4]), 0.5) # mean
vphi_l = c(p3[1]^2/3, p3[2]^2/3, (p3[3]^2*p5[3]+p4[3]^2*(1-p5[3]))-mphi_l[3]^2, (p3[4]^2*p5[4]+p4[4]^2*(1-p5[4]))-mphi_l[4]^2, 0) # variance
pm = rbind(mphi_l,vphi_l,p3,p4,p5)

# parameters of initial processes
init_l= c(1,2) # 2 cases
m0_l = c(100,3,1) # values of M_i

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
GenRandomPara = function(N,NumPara,case,ka=1) { ## DGPs for the working paper version in 2023
  para = array(0,dim=c(N,NumPara))
  ## initial values
  para[,1] = ka*0 # y1_{i} = 0 for all individuals
  ## AR(1) coefficients
  if (case==1) {
    mphi = 0.4; a = 0.3 # 0.4 + uniform (-0.3, 0.3)
    para[,3] = mphi + runif(N, min=-a, max=a) 
  }
  if (case==2) { 
    mphi = 0.4; a = 0.5 # 0.4 + uniform (-0.5, 0.5)
    para[,3] = mphi + runif(N, min=-a, max=a) 
  }
  if (case==3) {
    lb= 0.2; ub= 0.8; pi= 0.3 # categorical
    para[,3] = matrix(t(matrix(c(rep(lb,pi*N),rep(ub,(1-pi)*N)),N/10,10)),N,1)
  }
  if (case==4) {
    mphi = 0.475 # homogeneous phi_i = phi_0 = 0.475
    para[,3] = rep(mphi,N)
  }
  if (case==5) {
    mphi = 0.620 # homogeneous phi_i = phi_0 = 0.620
    para[,3] = rep(mphi,N)
  }
  ## individual fixed effects
  para[,2] = para[,3]+rnorm(N,0,1) 
  para
}

### Generate errors and a panel of outcome variables
GenError = function(N,T,ve,GAUSSIAN,GARCH){
  u = array(0,dim=c(N,T))
  e = array(0,dim=c(N,T))
  h = array(0,dim=c(N,T))
  
  if (GAUSSIAN == TRUE & GARCH == FALSE){
    e[,2:T] = matrix(rnorm(N*(T-1)),N,T-1)
    h[,2:T] = matrix(sqrt(ve),N,T-1)
    u[,2:T] = h[,2:T]*e[,2:T]
  }
  
  if (GAUSSIAN == TRUE & GARCH == TRUE){
    e[,2:T] = matrix(rnorm(N*(T-1)),N,T-1)
    cvs = t(matrix(mcvs,2,N)) # coefficients of GARCH effects
    h[,1] = sqrt(ve) # stationary initial values for GARCH(1,1)
    for(t in 2:T){
      ### GARCH(1,1)
      h[,t] = sqrt(ve-ve*rowSums(cvs) + cvs[,1]*h[,t-1]^2 + cvs[,2]*u[,t-1]^2)
      u[,t] = h[,t]*e[,t]
    }
  }
  
  if (GAUSSIAN == FALSE & GARCH == FALSE){
    e[,2:T] = matrix((rchisq(N*(T-1),2,ncp=0)-2)/2,N,T-1)
    h[,2:T] = matrix(sqrt(ve),N,T-1)
    u[,2:T] = h[,2:T]*e[,2:T]
  }
  
  if (GAUSSIAN == FALSE & GARCH == TRUE){
    cvs = t(matrix(mcvs,2,N)) # coefficients of GARCH effects
    e[,2:T] = matrix((rchisq(N*(T-1),2,ncp=0)-2)/2,N,T-1)
    h[,1] = sqrt(ve) # stationary initial values for GARCH(1,1)
    for(t in 2:T){
      ### GARCH(1,1)
      h[,t] = sqrt(ve-ve*rowSums(cvs) + cvs[,1]*h[,t-1]^2 + cvs[,2]*u[,t-1]^2)
      u[,t] = h[,t]*e[,t]
    }
  }
  
  return(u[,2:T])
}

### Generate Y
GenRandomY = function(N,Tmax,M0,para,ve,GAUSSIAN,GARCH){
  T = Tmax + M0 + 2
  Y = array(0,dim=c(N,T))
  Y[1:N,1] = para[1:N,1]
  GAUSSIAN = GAUSSIAN
  GARCH = GARCH
  errors = GenError(N,T,ve,GAUSSIAN,GARCH) # N*(T-1) matrix
  for(t in 1:(T-1)){
    ### Below: Y_{t+1} = alpha + beta Y_t + sigma 
    Y[1:N,(t+1)] = para[1:N,2] + para[1:N,3] * Y[1:N,t] + errors[1:N,t]
  }
  Y[,(M0+2+1):(M0+2+Tmax)]
}
##############################################################################

##############################################################################
#### FDAC and HetroGMM estimators of moments of AR(1) coefficients by Pesaran and Yang (2023)
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

MonteCarloSim = function(r,N,Tobs){
  panel0 = matrix(data_all[,r],Nmax,Tmax); 
  panel1 = panel0[1:N,1:Tobs]
  Moment_mv(panel1)
}

##############################################################################


##############################################################################
#### Exp1: Estimation of sigma_phi^2 in heterogeneous AR panels with uniformly distributed phi_i under the baseline case
##############################################################################
set.seed(123987654)
init = 1; b = 1; kappa = 2; GAUSSIAN = 1; GARCH = 0; m0 = 100; rep = 5000; reps = 1:rep
for (case in 1:2) {
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
      mc_results = mclapply(reps,MonteCarloSim,N=N,Tobs=Tobs,mc.preschedule = TRUE, mc.cores = numCores)

      tm = Tobs-1
      if (tm==3){
        mm_estimates = t(matrix(unlist(mc_results),4,rep))
        results_m1 = mm_estimates[,1]
        results_m1_sd = mm_estimates[,2]
      } # length = 4
      if (tm==4) {
        mm_estimates = t(matrix(unlist(mc_results),12,rep))
        results_m1 = mm_estimates[,1:2]
        results_m1_sd = mm_estimates[,3:4]
        results_m2 = mm_estimates[,5:6]
        results_m2_sd = mm_estimates[,7:8]
        results_var = mm_estimates[,9:10]
        results_var_sd = mm_estimates[,11:12]
      } # length = 4*3 = 12
      if (tm>4){
        mm_estimates = t(matrix(unlist(mc_results),16,rep))
        results_m1 = mm_estimates[,1:2]
        results_m1_sd = mm_estimates[,3:4]
        results_m2 = mm_estimates[,5:6]
        results_m2_sd = mm_estimates[,7:8]
        results_var = mm_estimates[,9:10]
        results_var_sd = mm_estimates[,11:12]
        results_m3 = mm_estimates[,13:14]
        results_m3_sd = mm_estimates[,15:16]
      }
      rm(mc_results);rm(list=ls(pattern="results*"));
      ###############################################################
      name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
      Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
      s = (Tid-1)*length(Nlist)+Nid;
      if (s<10) {
        name = paste(name,s,sep="_0")
      } else {
        name = paste(name,s,sep="_")
      }
      assign(name,mm_estimates);rm(name,mm_estimates)
      ###############################################################
    }
  }
  rm(data_all)
  fn = paste("exp_fdac_var_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
  save.image(file=fn)
  rm(list = ls(pattern = paste("^","mm_c",sep="")))
}
##############################################################################

##############################################################################
#### Exp2: Estimation of sigma_phi^2 in heterogeneous AR panels with uniformly distributed phi_i under different error processes
##############################################################################
init = 1; b = 1; kappa = 2; m0 = 100; rep = 5000; reps = 1:rep
for (g in 2:dim(ng)[1]) {
  GAUSSIAN = ng[g,1]; GARCH = ng[g,2]
  set.seed(123987654)
  for (case in 1:2) {
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
        mc_results = mclapply(reps,MonteCarloSim,N=N,Tobs=Tobs,mc.preschedule = TRUE, mc.cores = numCores)
        
        tm = Tobs-1
        if (tm==3){
          mm_estimates = t(matrix(unlist(mc_results),4,rep))
          results_m1 = mm_estimates[,1]
          results_m1_sd = mm_estimates[,2]
        } # length = 4
        if (tm==4) {
          mm_estimates = t(matrix(unlist(mc_results),12,rep))
          results_m1 = mm_estimates[,1:2]
          results_m1_sd = mm_estimates[,3:4]
          results_m2 = mm_estimates[,5:6]
          results_m2_sd = mm_estimates[,7:8]
          results_var = mm_estimates[,9:10]
          results_var_sd = mm_estimates[,11:12]
        } # length = 4*3 = 12
        if (tm>4){
          mm_estimates = t(matrix(unlist(mc_results),16,rep))
          results_m1 = mm_estimates[,1:2]
          results_m1_sd = mm_estimates[,3:4]
          results_m2 = mm_estimates[,5:6]
          results_m2_sd = mm_estimates[,7:8]
          results_var = mm_estimates[,9:10]
          results_var_sd = mm_estimates[,11:12]
          results_m3 = mm_estimates[,13:14]
          results_m3_sd = mm_estimates[,15:16]
        }
        rm(mc_results);rm(list=ls(pattern="results*"));
        ###############################################################
        name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
        Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
        s = (Tid-1)*length(Nlist)+Nid;
        if (s<10) {
          name = paste(name,s,sep="_0")
        } else {
          name = paste(name,s,sep="_")
        }
        assign(name,mm_estimates);rm(name,mm_estimates)
        ###############################################################
      }
    }
    rm(data_all)
    fn = paste("exp_fdac_var_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
    save.image(file=fn)
    rm(list = ls(pattern = paste("^","mm_c",sep="")))
  }
}
##############################################################################

