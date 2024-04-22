#### FDAC and HetroGMM estimators proposed in Pesaran and Yang (2024) "Hetro. AR panels"
rm(list=ls())

## Please change the path to the directory where the simulation outcomes are stores. 
# setwd('~/Downloads') 

## Please install the following packages. 
# list.of.packages <- c("openxlsx","latex2exp")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(openxlsx)
library(latex2exp)

#### Compare FDAC and HetroGMM estimators of sigma_phi^2 = Var(phi_i)

### Define Functions Calculate Statistics: (Bias, RMSE, T_test & Power, Coverage)
################################################################################
Bias = function(ests,truev) {
  ests = as.matrix(ests)
  truematrix = matrix(truev,nrow=nrow(ests),ncol=ncol(ests))
  colMeans(ests - truematrix)
}

RMSE = function(ests,truev) {
  ests = as.matrix(ests)
  truematrix = matrix(truev,nrow=nrow(ests),ncol=ncol(ests))
  sqrt(colMeans((ests - truematrix)^2))
}

T_test = function(truev,estsd,rep){
  ests = as.matrix(estsd[1:rep])
  sds = as.matrix(estsd[(rep+1):(2*rep)])
  truematrix = matrix(truev,nrow=nrow(ests),ncol=ncol(ests))
  colMeans(abs((ests-truematrix)/sds)>qnorm(p = 0.975))
}

Power = function(truev,estsd,pwgrid,rep){
  h1 = pwgrid
  powers = sapply(h1,function(x) T_test(x,estsd,rep))
  powers
}

Coverage = function(truev,estsd,rep){
  ests = as.matrix(estsd[1:rep])
  sds = as.matrix(estsd[(rep+1):(2*rep)])
  truem = matrix(truev,nrow=nrow(ests),ncol=ncol(ests))
  cvr = abs(truem - ests) < qnorm(0.975,0,1)*sds # within the 95% CI excluding 2.5% from each tail
  mean(cvr)
}
##############################################################################

## Create a style for center-aligning the numbers
##############################################################################
center_style <- createStyle(halign = "center")
right <- createStyle(halign = "right")
left <- createStyle(halign = "left")
bbs <- createStyle(border = "bottom",borderStyle="thin")
lbs <- createStyle(border = "bottom",borderStyle="double")
tbs <- createStyle(border = "top",borderStyle="double")
wrap_text_style <- createStyle(wrapText = TRUE)

k=2;
# parameters of AR(1) coefficients
p3 = c(0.5, 0.5, 0.5, 0.5, 0, 0) 
p4 = c(0, 0, 0.8, 1, 0, 0)
p5 = c(0, 0, 0.85, 0.95, 0, 0)
mphi_l = c(0.4, 0.5, p3[3]*p5[3]+p4[3]*(1-p5[3]), p3[4]*p5[4]+p4[4]*(1-p5[4]), 0.4, 0.5) # mean
vphi_l = c(p3[1]^2/3, p3[2]^2/3, (p3[3]^2*p5[3]+p4[3]^2*(1-p5[3]))-mphi_l[3]^2, (p3[4]^2*p5[4]+p4[4]^2*(1-p5[4]))-mphi_l[4]^2, 0, 0) # variance
pm = rbind(mphi_l,vphi_l,p3,p4,p5)
##############################################################################

Nlist2 = c(100,1000,2500,5000); Tlist2 = c(5,6,10); 
NT2 = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])

## Create a new workbook
wb <- createWorkbook()


#### Table: frequency of negative estimates
##############################################################################
var_np_1 = matrix(,2*dim(NT2)[1],2) ## Frequency of estimates being not positive over 2,000 replications

for (cid in 1:2) {
  GAUSSIAN = 1; GARCH = 0;
  fn = paste("exp_fdac_var_c",cid,"_i1_m100_b1_kp2_",GAUSSIAN,GARCH,".RData",sep="")
  load(fn)
  for (k in 1:dim(NT2)[1]){
    Tobs = NT2[k,1]; N = NT2[k,2];
    name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
    Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
    s = (Tid-1)*length(Nlist)+Nid;
    if (s<10) {
      name = paste(name,s,sep="_0")
    } else {
      name = paste(name,s,sep="_")
    }
    mc_results1 = get(name);rm(name)
    for (j in 1:2) { ## FDAC and HetroGMM estimators of sigma_phi^2
      if (Tobs>4) { 
        ## Frequency of estimates being not positive over the first 2,000 replications
        results_var1 = mc_results1[1:2000,c(9,11)]
        var_np_1[(cid-1)*dim(NT2)[1]+k,j] = mean(results_var1[,j]<0.0001) 
        rm(results_var1)
      }
    }
  }
}
rm(list = ls(pattern = paste("^","mm_c",sep="")));rm(mc_results1)

sn = "Table S.3"
addWorksheet(wb, sn)
writeData(wb, sn, x = "Table S.3: Frequency of FDAC and HetroGMM estimators of sigma_phi^2 = Var(phi_i) being negative with uniformly distributed phi_i and Gaussian errors without GARCH effects", startCol = 1, startRow = 1, colNames = FALSE)
writeData(wb, sn, x = "sigma_phi^2 = 0.083 with |phi_i|<1", startCol = 4, startRow = 2, colNames = FALSE)
writeData(wb, sn, x = "sigma_phi^2 = 0.083 with phi_i in [-1 + epsilon, 1] (epsilon>0)", startCol = 7, startRow = 2, colNames = FALSE)

Nlist2 = c(100,1000,2500,5000); Tlist2 = c(5,6,10); 
NT2 = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])
nr2 = dim(NT2)[1]; rl2 = seq(1,nr2,by=1);
Nlist2 = rep(c('100','1,000',"2,500",'5,000'),3)
nc = 2
tab = matrix(,nr2,2+(nc+1)*2)
tab2 = matrix(,nr2,2+(nc+1)*2)

#### Put numbers into the respective cells
for (r2 in 1:nr2) {
  tab2[r2,1] = NT2[r2,1]
  tab2[r2,2] = Nlist2[r2]
  
  tab[r2,4] = var_np_1[r2,1]
  tab[r2,5] = var_np_1[r2,2]
  tab[r2,7] = var_np_1[r2+12,1]
  tab[r2,8] = var_np_1[r2+12,2]
}

#### Change the formats of each column
for (col in c(4,5,7,8)) {   # frequency (*100)
  v0 = tab[,col]; v1 = v0
  for (i in 1:length(v0)) {
    if (is.na(v0[i]) == 0) {
      v1[i] = format(round(v0[i]*100, 1), nsmall = 1)
    } else {
      v1[i] = ""
    }
  }
  tab2[,col] = v1; rm(v0,v1)
}

h = tab2
h1 = c("T","n",rep(c("","FDAC","HetroGMM"),2))
rv0 = matrix(,1,8)
h = rbind(h1,h[1:4,],rv0,h[5:8,],rv0,h[9:12,])
rownames(h) <- NULL
colnames(h) <- NULL
writeData(wb,sn, x = h, startCol = 1, startRow = 3,colNames = FALSE, rowNames = FALSE)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
mergeCells(wb, sheet = sn, cols = 4:5, rows = 2)
mergeCells(wb, sheet = sn, cols = 7:8, rows = 2)

addStyle(wb,sn,style = center_style, rows = 2:3,cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = right, rows = 4:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:5,7:8), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 3,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = lbs, rows = (nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T,stack=T)
rm(tab,tab2)
##############################################################################


# Load the first 2,000 simulation results but only report summary statistics for the positive estimates
################################################################################
fn1 = "exp_fdac_var_c1_i1_m100_b1_kp2_10.RData"
fn2 = "exp_fdac_var_c2_i1_m100_b1_kp2_10.RData"; m0 = 100; b = 0; kp = 2
en = 1; temp1 = tryCatch(load(fn1),error=function(e) print(en))
temp2 = tryCatch(load(fn2),error=function(e) print(en))
case_list = c(1,2);

Nlist2 = c(100,1000,2500,5000); Tlist2 = c(5,6,10); 
NT2 = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])

if ( temp1[1]!=1 & temp2[1]!=1) {
  
  var_bias = matrix(,2*dim(NT2)[1],2)
  var_rmse = matrix(,2*dim(NT2)[1],2)
  var_size = matrix(,2*dim(NT2)[1],2)
  var_pw = matrix(,2*dim(NT2)[1],2*501)
  
  for (cid in 1:2) {
    
    case = case_list[cid]
    mb = mphi_l[case]
    vb = vphi_l[case]
    #### Get the respective simulation results
    for (k in 1:dim(NT2)[1]){
      Tobs = NT2[k,1]; N = NT2[k,2];
      name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
      Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
      s = (Tid-1)*length(Nlist)+Nid;
      if (s<10) {
        name = paste(name,s,sep="_0")
      } else {
        name = paste(name,s,sep="_")
      }
      mc_results1 = get(name);rm(name)
      for (j in 1:2) { # FDAC and HetroGMM estimators
        if (Tobs>4) { 
          results_var = mc_results1[1:2000,c(9,11)]
          results_var_sd = mc_results1[1:2000,c(10,12)]
          vid = which(results_var[,j]>=0.0001)
          rep = length(vid)
          var_bias[(cid-1)*dim(NT2)[1]+k,j] = Bias(results_var[vid,j], vb)
          var_rmse[(cid-1)*dim(NT2)[1]+k,j] = RMSE(results_var[vid,j], vb)
          v1 = cbind(results_var[vid,j],results_var_sd[vid,j])
          pwgrid = as.matrix(sort(c(seq(0.001,0.5,by=0.001),vb))) 
          var_pw[(cid-1)*dim(NT2)[1]+k,(501*(j-1)+1):(501*(j-1)+501)] = Power(vb,v1,pwgrid,rep = rep) # length = no. alternative values
          var_size[(cid-1)*dim(NT2)[1]+k,j] = T_test(vb,v1,rep = rep); rm(v1,vid)
        }
      }
    }
  }
  
  #### Table
  #########################################################################
  sn = "Table S.4"
  addWorksheet(wb, sn)
  writeData(wb, sn, x = "Table S.4: Bias, RMSE, and size of FDAC and HetroGMM estimators of sigma_phi^2 = Var(phi_i) in a heterogeneous panel AR(1) model with uniformly distributed phi_i and Gaussian errors without GARCH effects", startCol = 1, startRow = 1, colNames = FALSE)
  
  Nlist2 = c(100,1000,2500,5000); Tlist2 = c(5,6,10); 
  NT2 = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])
  nr2 = dim(NT2)[1]; rl2 = seq(1,nr2,by=1);
  Nlist2 = rep(c('100','1,000',"2,500",'5,000'),3)
  nc = 2
  tab = matrix(,nr2,2+(nc+1)*3*2)
  tab2 = matrix(,nr2,2+(nc+1)*3*2)
  
  #### Put numbers into the respective cells
  for (r2 in 1:nr2) {
    tab2[r2,1] = NT2[r2,1]
    tab2[r2,2] = Nlist2[r2]
    tab[r2,(3+1):(3+nc)] = var_bias[r2,]
    tab[r2,(3+nc+1+1):(3+nc+1+nc)] = var_rmse[r2,]
    tab[r2,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = var_size[r2,]
    
    tab[r2,(6+3*nc+1):(6+3*nc+nc)] = var_bias[r2+(nr2)*1,]
    tab[r2,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = var_rmse[r2+(nr2)*1,]
    tab[r2,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = var_size[r2+(nr2)*1,]
  }
  
  #### Change the formats of each column
  for (col in c(4:5,7:8,13:14,15:17)) {   # bias, RMSE
    v0 = tab[,col]; v1 = v0
    for (i in 1:length(v0)) {
      if (is.na(v0[i]) == 0) {
        v1[i] = format(round(v0[i], 3), nsmall = 3)
      } else {
        v1[i] = ""
      }
    }
    tab2[,col] = v1; rm(v0,v1)
  }
  for (col in c(10:11,19:20)) {   # size (*100)
    v0 = tab[,col]; v1 = v0
    for (i in 1:length(v0)) {
      if (is.na(v0[i]) == 0) {
        v1[i] = format(round(v0[i]*100, 1), nsmall = 1)
      } else {
        v1[i] = ""
      }
    }
    tab2[,col] = v1; rm(v0,v1)
  }
  
  h = tab2
  h1 = c("","","","Bias","","","RMSE","","","Size (*100)","","","Bias","","","RMSE","","","Size (*100)","")
  h2 = c("T","n","",rep(c('FDAC',"HetroGMM",""),5),"FDAC","HetroGMM")
  h3 = c(rep("",3),"sigma_phi^2 = 0.083 with |phi_i|<1",rep("",7+1),"sigma_phi^2 = 0.083 with phi_i in [-1 + epsilon, 1] for some epsilon > 0",rep("",7))
  rv0 = matrix(,1,20)
  h = rbind(h1,h2,h3,h[2:4,],rv0,h[6:8,],rv0,h[10:12,])
  rownames(h) <- NULL
  colnames(h) <- NULL
  writeData(wb,sn, x = h, startCol = 1, startRow = 2,colNames = FALSE, rowNames = FALSE)
  mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
  mergeCells(wb, sheet = sn, cols = 4:5, rows = 2)
  mergeCells(wb, sheet = sn, cols = 7:8, rows = 2)
  mergeCells(wb, sheet = sn, cols = 10:11, rows = 2)
  mergeCells(wb, sheet = sn, cols = 13:14, rows = 2)
  mergeCells(wb, sheet = sn, cols = 16:17, rows = 2)
  mergeCells(wb, sheet = sn, cols = 19:20, rows = 2)
  mergeCells(wb, sheet = sn, cols = 4:11, rows = 4)
  mergeCells(wb, sheet = sn, cols = 13:20, rows = 4)

  
  addStyle(wb,sn,style = center_style, rows = 2:4,cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = right, rows = 5:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:5,7:8,10:11,13:14,16:17,19:20), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 3,cols = 1:(ncol(h)), gridExpand = T,stack=T)  
  addStyle(wb,sn,style = bbs, rows = 4,cols = c(4:11,13:20), gridExpand = T,stack=T)
  addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = lbs, rows = (nrow(h)+1),cols = 1:(ncol(h)), gridExpand = T,stack=T)
  #########################################################################
  rm(tab,tab2)
  
  #### Figures 
  #########################################################################
  for (id in c(2)) {
    case = case_list[id]
    # Var(phi_i)
    #####################################################################################
    powermmv = t(var_pw[(1+(id-1)*dim(NT2)[1]):(id*dim(NT2)[1]),1:501]) # FDAC
    powergmmv = t(var_pw[(1+(id-1)*dim(NT2)[1]):(id*dim(NT2)[1]),502:1002]) # HetroGMM
    if (case==1) {
      name = "pw_fdac_hetrogmm_u1_a_var.png"
    }
    if (case==2) {
      name = "Figure S.2.png"
    }
    truev = vphi_l[case]
    vb = vphi_l[case]
    mcm = as.matrix(sort(c(seq(0.001,0.5,by=0.001),vb))) # values of alternatives (v0-0.5,v0+0.5) length = 101
    d1=10; d2 = 501
    mx = mcm[d1:d2]
    xmin = 0
    xmax = mx[length(mx)]
    ymin = 0
    ymax = 1
    
    png(name, units="in", width=30, height=38, res=50)
    par(mai=c(2,1.5,1.5,1),xpd=TRUE, mfrow=c(3,2),oma=c(6,3,3,2)) # (b,l,t,r)
    
    #### T=5
    plot(mx, powermmv[d1:d2,2], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="",main=expression(paste("FDAC estimator of ",sigma[phi]^2," (T=5)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powermmv[d1:d2,3], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powermmv[d1:d2,4], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.05), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(sigma[phi]^2), side = 1, line = 10, at = xmax , cex = 3) # label of x-axis
    mtext(expression(sigma[paste(phi,",0")]^2), side = 1, line = 10, at = truev , cex = 3) # label of x-axis    
    
    plot(mx, powergmmv[d1:d2,2], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="",main=expression(paste("HetroGMM estimator of ",sigma[phi]^2," (T=5)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powergmmv[d1:d2,3], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powergmmv[d1:d2,4], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.05), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(sigma[phi]^2), side = 1, line = 10, at = xmax , cex = 3) # label of x-axis
    mtext(expression(sigma[paste(phi,",0")]^2), side = 1, line = 10, at = truev , cex = 3) # label of x-axis  
    
    ### T=6
    plot(mx, powermmv[d1:d2,6], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("FDAC (T=6)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powermmv[d1:d2,7], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powermmv[d1:d2,8], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.05), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(sigma[phi]^2), side = 1, line = 10, at = xmax , cex = 3) # label of x-axis
    mtext(expression(sigma[paste(phi,",0")]^2), side = 1, line = 10, at = truev , cex = 3) # label of x-axis    
    
    plot(mx, powergmmv[d1:d2,6], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("HetroGMM (T=6)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powergmmv[d1:d2,7], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powergmmv[d1:d2,8], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.05), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(sigma[phi]^2), side = 1, line = 10, at = xmax , cex = 3) # label of x-axis
    mtext(expression(sigma[paste(phi,",0")]^2), side = 1, line = 10, at = truev , cex = 3) # label of x-axis    
    
    ### T=10
    plot(mx, powermmv[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("FDAC (T=10)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powermmv[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powermmv[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.05), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(sigma[phi]^2), side = 1, line = 10, at = xmax , cex = 3) # label of x-axis
    mtext(expression(sigma[paste(phi,",0")]^2), side = 1, line = 10, at = truev , cex = 3) # label of x-axis    
    
    plot(mx, powergmmv[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("HetroGMM (T=10)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powergmmv[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powergmmv[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.05), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(sigma[phi]^2), side = 1, line = 10, at = xmax , cex = 3) # label of x-axis
    mtext(expression(sigma[paste(phi,",0")]^2), side = 1, line = 10, at = truev , cex = 3) # label of x-axis    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottom',legend = c("n=1,000", "n=2,500", "n=5,000"), lty = c(4,2,1),col = c(1,2,4), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 5.5, seg.len=5, bty = 'n')
    
    dev.off()
    #####################################################################################
  }
  #########################################################################
  rm(list = ls(pattern = paste("^","mm_c",sep="")))
} 
################################################################################
rm(var_bias,var_rmse,var_size,var_pw)



# Load simulation results
################################################################################
case = 2; init = 1; m0 = 100;
mb = mphi_l[case]; vb = vphi_l[case]
Nlist2 = c(1000,2500,5000); Tlist2 = c(5,6,10); 
NT2 = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])

var_bias = matrix(,4*dim(NT2)[1],2)
var_rmse = matrix(,4*dim(NT2)[1],2)
var_size = matrix(,4*dim(NT2)[1],2)
var_pw = matrix(,4*dim(NT2)[1],2*501)
var_np = matrix(,4*dim(NT2)[1],2) ## Frequency of estimates being not positive over 2,000 replications

ng2 = matrix(c(1,0,1,0,0,0,1,1),4,2)

for (eid in 1:dim(ng2)[1]) {
  GAUSSIAN = ng2[eid,1]; GARCH = ng2[eid,2]
  fn = paste("exp_fdac_var_c",case,"_i",init,"_m",m0,"_b1_kp2_",GAUSSIAN,GARCH,".RData",sep="")
  en = -99999; temp = tryCatch(load(fn),error=function(e) print(en)); 
  
  if ( temp[1] != en ) {
    #### Get the respective simulation results
    for (k in 1:dim(NT2)[1]){
      Tobs = NT2[k,1]; N = NT2[k,2];
      name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
      Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
      s = (Tid-1)*length(Nlist)+Nid;
      if (s<10) {
        name = paste(name,s,sep="_0")
      } else {
        name = paste(name,s,sep="_")
      }
      mc_results1 = get(name);rm(name)
      for (j in 1:2) { # FDAC and HetroGMM estimators
        if (Tobs > 4) {
          results_var = mc_results1[1:2000,c(9,11)]
          results_var_sd = mc_results1[1:2000,c(10,12)]
          vid = which(results_var[,j]>=0.0001)
          rep = length(vid)
          var_np[(eid-1)*dim(NT2)[1]+k,j] =  1 - rep/2000
          var_bias[(eid-1)*dim(NT2)[1]+k,j] = Bias(results_var[vid,j], vb)
          var_rmse[(eid-1)*dim(NT2)[1]+k,j] = RMSE(results_var[vid,j], vb)
          v1 = cbind(results_var[vid,j],results_var_sd[vid,j])
          pwgrid = as.matrix(sort(c(seq(0.001,0.5,by=0.001),vb))) 
          var_pw[(eid-1)*dim(NT2)[1]+k,(501*(j-1)+1):(501*(j-1)+501)] = Power(vb,v1,pwgrid,rep = rep) # length = no. alternative values
          var_size[(eid-1)*dim(NT2)[1]+k,j] = T_test(vb,v1,rep = rep); rm(v1,vid)
          
        }
      }
    }
    rm(list = ls(pattern = paste("^","mm_c",sep="")))
  }
}

sn = "Table S.5"
addWorksheet(wb, sn)
writeData(wb, sn, x = "Table S.5: Bias, RMSE, and size of the FDAC estimator of sigma_phi^2 = 0.083 in a heterogeneous panel AR(1) model with uniformly distributed phi_i in [-1 + epsilon, 1] for some epsilon > 0 and different error processes", startCol = 1, startRow = 1, colNames = FALSE)

Nlist = c(1000,2500,5000); Tlist = c(5,6,10)
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
nr = dim(NT)[1]
nr2 = 2*nr; rl2 = seq(1,nr2,by=1);
Tlist2 = c("5","6","10"); Nlist2 = rep(c('1,000','2,500','5,000'),5)
nc = 2
tab = matrix(,nr2,2+(nc+1)*3)
tab2 = matrix(,nr2,2+(nc+1)*3)

#### Put numbers into the respective cells
tab = matrix(,nr2+1,2+(nc+1)*3*2)
tab2 = matrix(,nr2+1,2+(nc+1)*3*2)

#### Put numbers into the respective cells
for (r2 in 1:nr) {
  tab2[r2,1] = NT[r2,1]
  tab2[r2,2] = Nlist2[r2]
  tab[r2,(3+1):(3+nc)] = var_bias[r2,]
  tab[r2,(3+nc+1+1):(3+nc+1+nc)] = var_rmse[r2,]
  tab[r2,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = var_size[r2,]
  
  tab[r2,(6+3*nc+1):(6+3*nc+nc)] = var_bias[r2+nr*1,]
  tab[r2,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = var_rmse[r2+nr*1,]
  tab[r2,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = var_size[r2+nr*1,]
  
  tab2[r2+nr,1] = NT[r2,1]
  tab2[r2+nr,2] = Nlist2[r2]
  tab[r2+nr,(3+1):(3+nc)] = var_bias[r2+nr*2,]
  tab[r2+nr,(3+nc+1+1):(3+nc+1+nc)] = var_rmse[r2+nr*2,]
  tab[r2+nr,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = var_size[r2+nr*2,]
  
  tab[r2+nr,(6+3*nc+1):(6+3*nc+nc)] = var_bias[r2+nr*3,]
  tab[r2+nr,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = var_rmse[r2+nr*3,]
  tab[r2+nr,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = var_size[r2+nr*3,]
}

#### Change the formats of each column
for (col in c(4:5,7:8,13:14,16:17)) {   # bias, RMSE
  v0 = tab[,col]; v1 = v0
  for (i in 1:length(v0)) {
    if (is.na(v0[i]) == 0) {
      v1[i] = format(round(v0[i], 3), nsmall = 3)
    } else {
      v1[i] = ""
    }
  }
  tab2[,col] = v1; rm(v0,v1)
}
for (col in c(10:11,19:20)) {   # size (*100)
  v0 = tab[,col]; v1 = v0
  for (i in 1:length(v0)) {
    if (is.na(v0[i]) == 0) {
      v1[i] = format(round(v0[i]*100, 1), nsmall = 1)
    } else {
      v1[i] = ""
    }
  }
  tab2[,col] = v1; rm(v0,v1)
}

h = tab2
h1 = c("","","","Bias","","","RMSE","","","Size (*100)","","","Bias","","","RMSE","","","Size (*100)","")
h2 = c("T","n","",rep(c('FDAC',"HetroGMM",""),5),"FDAC","HetroGMM")
h3 = c(rep("",3),"Gaussian errors without GARCH effects",rep("",8),"Non-Gaussian errors without GARCH effects",rep("",7))
h4 = c(rep("",3),"Gaussian errors with GARCH effects",rep("",8),"Gaussian errors with GARCH effects",rep("",7))

rv0 = matrix(,1,20)
h = rbind(h1,h2,h3,h[1:3,],rv0,h[4:6,],rv0,h[7:9,],rv0,h4,h[10:12,],rv0,h[13:15,],rv0,h[16:18,])
rownames(h) <- NULL
colnames(h) <- NULL
writeData(wb,sn, x = h, startCol = 1, startRow = 2,colNames = FALSE, rowNames = FALSE)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
mergeCells(wb, sheet = sn, cols = 4:5, rows = 2)
mergeCells(wb, sheet = sn, cols = 7:8, rows = 2)
mergeCells(wb, sheet = sn, cols = 10:11, rows = 2)
mergeCells(wb, sheet = sn, cols = 13:14, rows = 2)
mergeCells(wb, sheet = sn, cols = 16:17, rows = 2)
mergeCells(wb, sheet = sn, cols = 19:20, rows = 2)
mergeCells(wb, sheet = sn, cols = 4:11, rows = 4)
mergeCells(wb, sheet = sn, cols = 13:20, rows = 4)
mergeCells(wb, sheet = sn, cols = 4:11, rows = 17)
mergeCells(wb, sheet = sn, cols = 13:20, rows = 17)

addStyle(wb,sn,style = right, rows = 4:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = center_style, rows = c(2:4,17),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:7,9:12,14:17), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 3,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 4,cols = c(4:11,13:20), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 17,cols = c(4:11,13:20), gridExpand = T,stack=T)
addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = lbs, rows = (nrow(h)+1),cols = 1:(ncol(h)), gridExpand = T,stack=T)
rm(tab,tab2)

sn = "Table S.7"
addWorksheet(wb, sn)
writeData(wb, sn, x = "Table S.7: Frequency of FDAC and HetroGMM estimators of sigma_phi^2 = Var(phi_i) being negative with uniformly distributed phi_i in [-1 + epsilon, 1] for some epsilon > 0 under different error processes", startCol = 1, startRow = 1, colNames = FALSE)
writeData(wb, sn, x = "Without GARCH effects", startCol = 4, startRow = 2, colNames = FALSE)
writeData(wb, sn, x = "With GARCH effects", startCol = 10, startRow = 2, colNames = FALSE)

Nlist2 = c(1000,2500,5000); Tlist2 = c(5,6,10); 
NT2 = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])
nr2 = dim(NT2)[1]; rl2 = seq(1,nr2,by=1);
Nlist2 = rep(c('1,000',"2,500",'5,000'),3)
nc = 5
tab = matrix(,nr2,2+(nc+1)*2)
tab2 = matrix(,nr2,2+(nc+1)*2)

#### Put numbers into the respective cells
for (r2 in 1:nr2) {
  tab2[r2,1] = NT2[r2,1]
  tab2[r2,2] = Nlist2[r2]
  
  tab[r2,4] = var_np[r2,1]
  tab[r2,5] = var_np[r2,2]
  
  tab[r2,7] = var_np[r2+9,1]
  tab[r2,8] = var_np[r2+9,2]
  
  tab[r2,10] = var_np[r2+18,1]
  tab[r2,11] = var_np[r2+18,2]
  
  tab[r2,13] = var_np[r2+27,1]
  tab[r2,14] = var_np[r2+27,2]
}

#### Change the formats of each column
for (col in c(4,5,7,8,10,11,13,14)) {   # frequency (*100)
  v0 = tab[,col]; v1 = v0
  for (i in 1:length(v0)) {
    if (is.na(v0[i]) == 0) {
      v1[i] = format(round(v0[i]*100, 1), nsmall = 1)
    } else {
      v1[i] = ""
    }
  }
  tab2[,col] = v1; rm(v0,v1)
}

h = tab2
h1 = c(rep("",2),rep(c("","Gaussian","","","Non-Gaussian",""),2))
h2 = c("T","n",rep(c("","FDAC","HetroGMM"),4))
rv0 = matrix(,1,14)
h = rbind(h1,h2,h[1:3,],rv0,h[4:6,],rv0,h[7:9,])
rownames(h) <- NULL
colnames(h) <- NULL
writeData(wb,sn, x = h, startCol = 1, startRow = 3,colNames = FALSE, rowNames = FALSE)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
mergeCells(wb, sheet = sn, cols = 4:8, rows = 2)
mergeCells(wb, sheet = sn, cols = 10:14, rows = 2)
mergeCells(wb, sheet = sn, cols = 4:5, rows = 3)
mergeCells(wb, sheet = sn, cols = 7:8, rows = 3)
mergeCells(wb, sheet = sn, cols = 10:11, rows = 3)
mergeCells(wb, sheet = sn, cols = 13:14, rows = 3)

addStyle(wb,sn,style = center_style, rows = 2:4,cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = right, rows = 5:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:8,10:14), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 3,cols = c(4:5,7:8,10:11,13:14), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 4,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = lbs, rows = (nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T,stack=T)
rm(tab,tab2)
################################################################################



# Save the workbook to an Excel file
saveWorkbook(wb, file = "HetroAR_MC_FDAC_moments_exp_var.xlsx",overwrite = TRUE)
cat("The MC results have been written to the excel file HetroAR_MC_FDAC_moments_exp_var.xlsx.")

