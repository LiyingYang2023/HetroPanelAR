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

#### Compare FDAC and HetroGMM estimators of mu_phi = E(phi_i)

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
##############################################################################

## Create a new workbook
wb <- createWorkbook()

# parameters of AR(1) coefficients
p3 = c(0.5, 0.5, 0.5, 0.5, 0, 0) 
p4 = c(0, 0, 0.8, 1, 0, 0)
p5 = c(0, 0, 0.85, 0.95, 0, 0)
mphi_l = c(0.4, 0.5, p3[3]*p5[3]+p4[3]*(1-p5[3]), p3[4]*p5[4]+p4[4]*(1-p5[4]), 0.4, 0.5) # mean
vphi_l = c(p3[1]^2/3, p3[2]^2/3, (p3[3]^2*p5[3]+p4[3]^2*(1-p5[3]))-mphi_l[3]^2, (p3[4]^2*p5[4]+p4[4]^2*(1-p5[4]))-mphi_l[4]^2, 0, 0) # variance
pm = rbind(mphi_l,vphi_l,p3,p4,p5)


## Table 1: Compare FDAC and HetroGMM estimators of mu_phi = E(phi_i) with uniformly distributed phi_i and Figure 1
##############################################################################
# Load simulation results
fn1 = "exp_fdac_mm_c1_i1_m100_b1_kp2_10.RData"
fn2 = "exp_fdac_mm_c2_i1_m100_b1_kp2_10.RData"; m0 = 100; b = 0; kp = 2; rep = 2000
en = 1; temp1 = tryCatch(load(fn1),error=function(e) print(en))
temp2 = tryCatch(load(fn2),error=function(e) print(en))
case_list = c(1,2);
Nlist = c(100,1000,2500,5000); 
Tlist = c(4,5,6,10); 
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)


if ( temp1[1]!=1 & temp2[1]!=1) {
  m1_bias = matrix(,2*length(NT_l),2)
  m1_rmse = matrix(,2*length(NT_l),2)
  m1_size = matrix(,2*length(NT_l),2)
  m1_pw = matrix(,2*length(NT_l),2*501)

  for (cid in 1:2) {

    case = case_list[cid]
    mb = mphi_l[case]
    #### Get the respective simulation results
    for (k in 1:length(NT_l)){
      kid = NT_l[k]
      Tobs = NT[kid,1]; N = NT[kid,2];
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
          results_m1 = mc_results1[1:2000,c(1,3)]
          results_m1_sd = mc_results1[1:2000,c(2,4)]
          m1_bias[(cid-1)*length(NT_l)+k,j] = Bias(results_m1[,j], mb)
          m1_rmse[(cid-1)*length(NT_l)+k,j] = RMSE(results_m1[,j], mb)
          e1 = cbind(results_m1[,j],results_m1_sd[,j])
          pwgrid = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01)) 
          m1_pw[(cid-1)*length(NT_l)+k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,pwgrid,rep=2000) # length = no. alternative values
          m1_size[(cid-1)*length(NT_l)+k,j] = T_test(mb,e1,rep=2000); rm(e1)
      }
    }
  }

  #### Table of E(phi_i)
  #########################################################################
  sn = "Table S.1"
  addWorksheet(wb, sn)
  writeData(wb, sn, x = "Table S.1: Bias, RMSE, and size of FDAC and HetroGMM estimators of mu_phi = E(phi_i) in a heterogeneous panel AR(1) model with uniformly distributed phi_i and Gaussian errors without GARCH effects", startCol = 1, startRow = 1, colNames = FALSE)
  writeData(wb, sn, x = "mu_phi = 0.4 with |phi_i|<1", startCol = 4, startRow = 2, colNames = FALSE)
  writeData(wb, sn, x = "mu_phi = 0.5 with phi_i in [-1 + epsilon, 1] Keywords:} Heterogeneous dynamic panels, neglected
heterogeneity bias, short $T$ panels, earnings dynamicsfor some epsilon > 0", startCol = 13, startRow = 2, colNames = FALSE)

  Nlist = c(100,1000,5000); Tlist = c(4,5,6,10)
  NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
  nr2 = dim(NT)[1]; rl2 = seq(1,nr2,by=1);
  Tlist2 = c("4","5","6","10"); Nlist2 = rep(c('100','1,000','5,000'),5)
  nc = 2
  tab = matrix(,nr2+1,2+(nc+1)*3*2)
  tab2 = matrix(,nr2+1,2+(nc+1)*3*2)

  #### Put numbers into the respective cells
  for (r2 in 1:nr2) {
    tab2[r2,1] = NT[r2,1]
    tab2[r2,2] = Nlist2[r2]
    tab[r2,(3+1):(3+nc)] = m1_bias[r2,]
    tab[r2,(3+nc+1+1):(3+nc+1+nc)] = m1_rmse[r2,]
    tab[r2,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = m1_size[r2,]

    tab[r2,(6+3*nc+1):(6+3*nc+nc)] = m1_bias[r2+nr2*1,]
    tab[r2,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = m1_rmse[r2+nr2*1,]
    tab[r2,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = m1_size[r2+nr2*1,]
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
  rv0 = matrix(,1,20)
  h = rbind(h1,h2,h[1:3,],rv0,h[4:6,],rv0,h[7:9,],rv0,h[10:12,])
  rownames(h) <- NULL
  colnames(h) <- NULL
  writeData(wb,sn, x = h, startCol = 1, startRow = 3,colNames = FALSE, rowNames = FALSE)
  mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
  mergeCells(wb, sheet = sn, cols = 4:11, rows = 2)
  mergeCells(wb, sheet = sn, cols = 13:20, rows = 2)
  mergeCells(wb, sheet = sn, cols = 4:5, rows = 3)
  mergeCells(wb, sheet = sn, cols = 7:8, rows = 3)
  mergeCells(wb, sheet = sn, cols = 10:11, rows = 3)
  mergeCells(wb, sheet = sn, cols = 13:14, rows = 3)
  mergeCells(wb, sheet = sn, cols = 16:17, rows = 3)
  mergeCells(wb, sheet = sn, cols = 19:20, rows = 3)

  addStyle(wb,sn,style = center_style, rows = 2:4,cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = right, rows = 5:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:11,13:20), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 3,cols = c(4:5,7:8,10:11,13:14,16:17,19:20), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 4,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = lbs, rows = (nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T,stack=T)
  #########################################################################
  rm(tab,tab2)


  #### Figures of mu_phi = E(phi_i) 
  #########################################################################
  for (id in 2) {
    case = case_list[id]
    # E(phi_i)
    #####################################################################################
    powermm1 = t(m1_pw[(1+(id-1)*length(NT_l)):(id*length(NT_l)),1:501]) # FDAC
    powergmm1 = t(m1_pw[(1+(id-1)*length(NT_l)):(id*length(NT_l)),502:1002]) # HetroGMM
    if (case==1) {
      name = "pw_fdac_hetrogmm_u1_a_m1.png"
    }
    if (case==2) {
      name = "Figure S.1.png"
    }
    #### Define x axis by alternative values
    d1 = 151; d2 = 351
    truev = mphi_l[case]
    mb = mphi_l[case]
    mcm = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01)) 
    mx = mcm[d1:d2]
    xmin = mx[1]
    xmax = mx[length(mx)]
    ymin = 0
    ymax = 1

    png(name, units="in", width=40, height=38, res=50)
    par(mai=c(2,1.5,1.5,1),xpd=TRUE, mfrow=c(3,2),oma=c(6,3,3,2)) # (b,l,t,r)

    #### T=4
    plot(mx, powermm1[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("FDAC estimator of ",mu[phi]," (T=4)")),cex.lab=4, axes = F,cex.main=5.5, cex.sub=4)
    lines(mx, powermm1[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powermm1[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.1), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(mu[phi]), side = 1, line = 9, at = xmax , cex = 3.5) # label of x-axis
    mtext(expression(mu[paste(phi,",0")]), side = 1, line = 9, at = truev , cex = 3.5) # label of x-axis    


    plot(mx, powergmm1[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="",main=expression(paste("HetroGMM estimator of ",mu[phi]," (T=4)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powergmm1[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powergmm1[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.1), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(mu[phi]), side = 1, line = 9, at = xmax , cex = 3.5) # label of x-axis
    mtext(expression(mu[paste(phi,",0")]), side = 1, line = 9, at = truev , cex = 3.5) # label of x-axis    
    
    ### T=6
    plot(mx, powermm1[d1:d2,7], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("FDAC (T=6)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powermm1[d1:d2,8], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powermm1[d1:d2,9], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.1), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(mu[phi]), side = 1, line = 9, at = xmax , cex = 3.5) # label of x-axis
    mtext(expression(mu[paste(phi,",0")]), side = 1, line = 9, at = truev , cex = 3.5) # label of x-axis    

    plot(mx, powergmm1[d1:d2,7], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("HetroGMM (T=6)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powergmm1[d1:d2,8], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powergmm1[d1:d2,9], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.1), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(mu[phi]), side = 1, line = 9, at = xmax , cex = 3.5) # label of x-axis
    mtext(expression(mu[paste(phi,",0")]), side = 1, line = 9, at = truev , cex = 3.5) # label of x-axis    

    ### T=10
    plot(mx, powermm1[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("FDAC (T=10)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powermm1[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powermm1[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.1), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(mu[phi]), side = 1, line = 9, at = xmax , cex = 3.5) # label of x-axis
    mtext(expression(mu[paste(phi,",0")]), side = 1, line = 9, at = truev , cex = 3.5) # label of x-axis    

    plot(mx, powergmm1[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="",main=expression(paste("HetroGMM (T=10)")),cex.lab=4,  axes = F, cex.main=5.5, cex.sub=4)
    lines(mx, powergmm1[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    lines(mx, powergmm1[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    axis(side=1, at=seq(xmin, xmax, by=0.1), mgp=c(2, 3.5, 0), cex.axis=3,labels = T) # x-axis
    axis(side=2, at=seq(0, 1, by=0.1), mgp=c(4, 2, 0), cex.axis=3) # y-axis
    arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
    arrows(1,-0.01, 1, 1, length = 0,lty =3,lwd=3) # denote mu_phi = 1
    arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
    mtext("power", side = 2, line = 6, at = 1 , cex = 3) # label of y-axis
    mtext("5%", side = 1, line = -6, at = xmax , cex = 3) 
    mtext(expression(mu[phi]), side = 1, line = 9, at = xmax , cex = 3.5) # label of x-axis
    mtext(expression(mu[paste(phi,",0")]), side = 1, line = 9, at = truev , cex = 3.5) # label of x-axis    

    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottom',legend = c("n=100", "n=1,000", "n=5,000"), lty = c(4,2,1),col = c(1,2,4), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 5.5, seg.len=5, bty = 'n')

    dev.off()
    #####################################################################################
  }
  #########################################################################
  rm(list = ls(pattern = paste("^","mm_c",sep="")))
} 
##############################################################################

## Table 2: Compare FDAC and HetroGMM estimators of mu_phi = E(phi_i) with categorically distributed phi_i
##############################################################################
# Load simulation results
fn1 = "exp_fdac_mm_c3_i1_m100_b1_kp2_10.RData"
fn2 = "exp_fdac_mm_c4_i1_m100_b1_kp2_10.RData"; m0 = 100; b = 0; kp = 2; rep = 2000
en = 1; temp1 = tryCatch(load(fn1),error=function(e) print(en))
temp2 = tryCatch(load(fn2),error=function(e) print(en))
case_list = c(3,4);
Nlist = c(100,1000,2500,5000); 
Tlist = c(4,5,6,10); 
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)


if ( temp1[1]!=1 & temp2[1]!=1) {
  m1_bias = matrix(,2*length(NT_l),2)
  m1_rmse = matrix(,2*length(NT_l),2)
  m1_size = matrix(,2*length(NT_l),2)
  m1_pw = matrix(,2*length(NT_l),2*501)
  
  for (cid in 1:2) {
    
    case = case_list[cid]
    mb = mphi_l[case]
    #### Get the respective simulation results
    for (k in 1:length(NT_l)){
      kid = NT_l[k]
      Tobs = NT[kid,1]; N = NT[kid,2];
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
        results_m1 = mc_results1[1:2000,c(1,3)]
        results_m1_sd = mc_results1[1:2000,c(2,4)]
        m1_bias[(cid-1)*length(NT_l)+k,j] = Bias(results_m1[,j], mb)
        m1_rmse[(cid-1)*length(NT_l)+k,j] = RMSE(results_m1[,j], mb)
        e1 = cbind(results_m1[,j],results_m1_sd[,j])
        pwgrid = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01)) 
        m1_pw[(cid-1)*length(NT_l)+k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,pwgrid,rep=2000) # length = no. alternative values
        m1_size[(cid-1)*length(NT_l)+k,j] = T_test(mb,e1,rep=2000); rm(e1)
      }
    }
  }
  
  #### Table of E(phi_i) 
  #########################################################################
  sn = "Table S.2"
  addWorksheet(wb, sn)
  writeData(wb, sn, x = "Table S.2: Bias, RMSE, and size of FDAC and HetroGMM estimators of mu_phi = E(phi_i) in a heterogeneous panel AR(1) model with categorically distributed phi_i and Gaussian errors without GARCH effects", startCol = 1, startRow = 1, colNames = FALSE)
  writeData(wb, sn, x = "mu_phi = 0.525 with |phi_i|<1", startCol = 4, startRow = 2, colNames = FALSE)
  writeData(wb, sn, x = "mu_phi = 0.545 with phi_i in [-1 + epsilon, 1] for some epsilon > 0", startCol = 13, startRow = 2, colNames = FALSE)
  
  Nlist = c(100,1000,5000); Tlist = c(4,5,6,10)
  NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
  nr2 = dim(NT)[1]; rl2 = seq(1,nr2,by=1);
  Tlist2 = c("4","5","6","10"); Nlist2 = rep(c('100','1,000','5,000'),5)
  nc = 2
  tab = matrix(,nr2+1,2+(nc+1)*3*2)
  tab2 = matrix(,nr2+1,2+(nc+1)*3*2)
  
  #### Put numbers into the respective cells
  for (r2 in 1:nr2) {
    tab2[r2,1] = NT[r2,1]
    tab2[r2,2] = Nlist2[r2]
    tab[r2,(3+1):(3+nc)] = m1_bias[r2,]
    tab[r2,(3+nc+1+1):(3+nc+1+nc)] = m1_rmse[r2,]
    tab[r2,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = m1_size[r2,]
    
    tab[r2,(6+3*nc+1):(6+3*nc+nc)] = m1_bias[r2+nr2*1,]
    tab[r2,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = m1_rmse[r2+nr2*1,]
    tab[r2,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = m1_size[r2+nr2*1,]
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
  rv0 = matrix(,1,20)
  h = rbind(h1,h2,h[1:3,],rv0,h[4:6,],rv0,h[7:9,],rv0,h[10:12,])
  rownames(h) <- NULL
  colnames(h) <- NULL
  writeData(wb,sn, x = h, startCol = 1, startRow = 3,colNames = FALSE, rowNames = FALSE)
  mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
  mergeCells(wb, sheet = sn, cols = 4:11, rows = 2)
  mergeCells(wb, sheet = sn, cols = 13:20, rows = 2)
  mergeCells(wb, sheet = sn, cols = 4:5, rows = 3)
  mergeCells(wb, sheet = sn, cols = 7:8, rows = 3)
  mergeCells(wb, sheet = sn, cols = 10:11, rows = 3)
  mergeCells(wb, sheet = sn, cols = 13:14, rows = 3)
  mergeCells(wb, sheet = sn, cols = 16:17, rows = 3)
  mergeCells(wb, sheet = sn, cols = 19:20, rows = 3)
  
  addStyle(wb,sn,style = center_style, rows = 2:4,cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = right, rows = 5:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:11,13:20), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 3,cols = c(4:5,7:8,10:11,13:14,16:17,19:20), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 4,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = lbs, rows = (nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T,stack=T)
  #########################################################################
  rm(tab,tab2)
  rm(list = ls(pattern = paste("^","mm_c",sep="")))
} 
##############################################################################

## Tables 3: Compare the FDAC estimator of mean with uniformly distributed phi_i under different error processes
##############################################################################
# Load simulation results
case = 2; init = 1; m0 = 100; b = 1; kp = 2
mb = mphi_l[case]; vb = vphi_l[case]
Nlist = c(100,1000,2500,5000); Tlist = c(4,5,6,10)
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)

m1_bias = matrix(,4*length(NT_l),2)
m1_rmse = matrix(,4*length(NT_l),2)
m1_size = matrix(,4*length(NT_l),2)
m1_pw = matrix(,4*length(NT_l),2*501)

ng2 = matrix(c(1,0,1,0,0,0,1,1),4,2)

for (eid in 1:dim(ng2)[1]) {
  GAUSSIAN = ng2[eid,1]; GARCH = ng2[eid,2]
  fn = paste("exp_fdac_mm_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
  en = -99999; temp = tryCatch(load(fn),error=function(e) print(en)); rep = 2000

  if ( temp[1] != en ) {
    #### Get the respective simulation results
    for (k in 1:length(NT_l)){
      kid = NT_l[k]
      Tobs = NT[kid,1]; N = NT[kid,2];
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
        results_m1 = mc_results1[1:2000,c(1,3)]
        results_m1_sd = mc_results1[1:2000,c(2,4)]
        m1_bias[(eid-1)*length(NT_l)+k,j] = Bias(results_m1[,j], mb)
        m1_rmse[(eid-1)*length(NT_l)+k,j] = RMSE(results_m1[,j], mb)
        e1 = cbind(results_m1[,j],results_m1_sd[,j])
        pwgrid = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01)) 
        m1_pw[(eid-1)*length(NT_l)+k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,pwgrid,rep=2000) # length = no. alternative values
        m1_size[(eid-1)*length(NT_l)+k,j] = T_test(mb,e1,rep=2000); rm(e1)
      }
    }
    rm(list = ls(pattern = paste("^","mm_c",sep="")))
  }
}

sn = "Table S.4"
addWorksheet(wb, sn)
writeData(wb, sn, x = "Table S.4: Bias, RMSE, and size of the FDAC estimator of mu_phi = E(phi_i) = 0.5 in a heterogeneous panel AR(1) model with uniformly distributed phi_i in [-1 + epsilon, 1] for some epsilon > 0 and different error processes", startCol = 1, startRow = 1, colNames = FALSE)

Nlist = c(100,1000,5000); Tlist = c(4,5,6,10)
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
nr = dim(NT)[1]
nr2 = 2*nr; rl2 = seq(1,nr2,by=1);
Tlist2 = c("4","5","6","10"); Nlist2 = rep(c('100','1,000','5,000'),5)
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
  tab[r2,(3+1):(3+nc)] = m1_bias[r2,]
  tab[r2,(3+nc+1+1):(3+nc+1+nc)] = m1_rmse[r2,]
  tab[r2,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = m1_size[r2,]
  
  tab[r2,(6+3*nc+1):(6+3*nc+nc)] = m1_bias[r2+nr*1,]
  tab[r2,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = m1_rmse[r2+nr*1,]
  tab[r2,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = m1_size[r2+nr*1,]
  
  tab2[r2+nr,1] = NT[r2,1]
  tab2[r2+nr,2] = Nlist2[r2]
  tab[r2+nr,(3+1):(3+nc)] = m1_bias[r2+nr*2,]
  tab[r2+nr,(3+nc+1+1):(3+nc+1+nc)] = m1_rmse[r2+nr*2,]
  tab[r2+nr,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = m1_size[r2+nr*2,]
  
  tab[r2+nr,(6+3*nc+1):(6+3*nc+nc)] = m1_bias[r2+nr*3,]
  tab[r2+nr,(6+3*nc+nc+1+1):(6+3*nc+nc+1+nc)] = m1_rmse[r2+nr*3,]
  tab[r2+nr,(6+3*nc+nc+1+nc+1+1):(6+3*nc+nc+1+nc+1+nc)] = m1_size[r2+nr*3,]
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
h4 = c(rep("",3),"Gaussian errors with GARCH effects",rep("",8),"Non-Gaussian errors with GARCH effects",rep("",7))

rv0 = matrix(,1,20)
h = rbind(h1,h2,h3,h[1:3,],rv0,h[4:6,],rv0,h[7:9,],rv0,h[10:12,],rv0,h4,h[13:15,],rv0,h[16:18,],rv0,h[19:21,],rv0,h[22:24,])
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
mergeCells(wb, sheet = sn, cols = 4:11, rows = 21)
mergeCells(wb, sheet = sn, cols = 13:20, rows = 21)

addStyle(wb,sn,style = right, rows = 4:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = center_style, rows = c(2:4,21),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:7,9:12,14:17), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 3,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 4,cols = c(4:11,13:20), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 21,cols = c(4:11,13:20), gridExpand = T,stack=T)
addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = lbs, rows = (nrow(h)+1),cols = 1:(ncol(h)), gridExpand = T,stack=T)
rm(tab,tab2)
##############################################################################

## Table 4: Compare the FDAC estimator of mu_phi = E(phi_i) under different initial processes
##############################################################################
# Load simulation results
sn_list = c("Table S.12", "Table S.13")
tabn_list = c("Table S.12: Bias, RMSE, and size of the FDAC estimator of mu_phi = E(phi_i) in a heterogeneous panel AR(1) model with uniformly distributed phi_i, Gaussian errors without GARCH effects, and different initializations",
              "Table S.13: Bias, RMSE, and size of the FDAC estimator of mu_phi = E(phi_i) in a heterogeneous panel AR(1) model with categorically distributed phi_i, Gaussian errors without GARCH effects, and different initializations")


Nlist = c(100,1000,2500,5000); Tlist = c(4,5,6,10)
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
nc = 3
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)
m1_bias = matrix(,4*length(NT_l),nc)
m1_rmse = matrix(,4*length(NT_l),nc)
m1_size = matrix(,4*length(NT_l),nc)
m1_pw = matrix(,4*length(NT_l),nc*501)

for (case in 1:4) {
  mb = mphi_l[case]; vb = vphi_l[case]
  
  for (j in 1:nc) {
    init = 1; b = 1; kappa = 2;
    m0_l = c(100,3,1)
    m0 = m0_l[j]
    GAUSSIAN = 1; GARCH = 0
    fn = paste("exp_fdac_mm_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
    en = -99999; temp = tryCatch(load(fn),error=function(e) print(en))
    rep = 2000
    
    Nlist = c(100,1000,2500,5000); Tlist = c(4,5,6,10)
    NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
    NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)
    nc = 3
    
    if ( temp[1] != en ) {
      #### Get the respective simulation results
      for (k in 1:length(NT_l)){
        kid = NT_l[k]
        Tobs = NT[kid,1]; N = NT[kid,2];
        name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
        Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
        s = (Tid-1)*length(Nlist)+Nid;
        if (s<10) {
          name = paste(name,s,sep="_0")
        } else {
          name = paste(name,s,sep="_")
        }
        mc_results1 = get(name);rm(name)
        
        if (Tobs>3) { # Mean
          results_m1 = mc_results1[1:2000,c(1,3)]
          results_m1_sd = mc_results1[1:2000,c(2,4)]
          m1_bias[k+(case-1)*12,j] = Bias(results_m1[,1], mb)
          m1_rmse[k+(case-1)*12,j] = RMSE(results_m1[,1], mb)
          e1 = cbind(results_m1[,1],results_m1_sd[,1])
          pwgrid = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01))
          m1_pw[k+(case-1)*12,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,pwgrid,rep) # length = no. alternative values
          m1_size[k+(case-1)*12,j] = T_test(mb,e1,rep); rm(e1)
        }
      }
      rm(list = ls(pattern = paste("^","mm_c",sep="")))
    }
  }
}
  

for (d in 1:2) {
  sn = sn_list[d]
  addWorksheet(wb, sn)
  writeData(wb, sn, x = tabn_list[d], startCol = 1, startRow = 1, colNames = FALSE)
  writeData(wb, sn, x = "Bias", startCol = 4, startRow = 2, colNames = FALSE)
  writeData(wb, sn, x = "RMSE", startCol = 8, startRow = 2, colNames = FALSE)
  writeData(wb, sn, x = "Size (*100)", startCol = 12, startRow = 2, colNames = FALSE)
  
  Nlist = c(100,1000,5000); Tlist = c(4,5,6,10)
  NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
  nr2 = dim(NT)[1]; rl2 = seq(1,nr2,by=1);
  Tlist2 = c("4","5","6","10"); Nlist2 = rep(c('100','1,000','5,000'),4)
  tab = matrix(,nr2*2,2+(nc+1)*3)
  tab2 = matrix(,nr2*2,2+(nc+1)*3)
  
  #### Put numbers into the respective cells
  for (r2 in 1:nr2) {
    tab2[r2,1] = NT[r2,1]
    tab2[r2,2] = Nlist2[r2]
    
    tab2[r2+12,1] = NT[r2,1]
    tab2[r2+12,2] = Nlist2[r2]
  }
  tab[1:12,(3+1):(3+nc)] = m1_bias[(1:12)+(d-1)*24,]
  tab[13:24,(3+1):(3+nc)] = m1_bias[(13:24)+(d-1)*24,]
  
  tab[1:12,(3+nc+1+1):(3+nc+1+nc)] = m1_rmse[(1:12)+(d-1)*24,]
  tab[13:24,(3+nc+1+1):(3+nc+1+nc)] = m1_rmse[(13:24)+(d-1)*24,]
  
  tab[1:12,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = m1_size[(1:12)+(d-1)*24,]
  tab[13:24,(3+nc+1+nc+1+1):(3+nc+1+nc+1+nc)] = m1_size[(13:24)+(d-1)*24,]
  
  #### Change the formats of each column
  for (col in c(4:6, 8:10)) {   # bias, RMSE
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
  for (col in c(12:14)) {   # size (*100)
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
  h2 = c("T","n/M_i",rep(c('','100',"3","1"),3))
  if (d == 1) {
    h3 = c("mu_phi = 0.4 with |phi_i| < 1",rep("",13))
    h4 = c("mu_phi = 0.5 with phi_i in [-1 + epsilon, 1] for some epsilon > 0",rep("",13))
  }
  if (d == 2) {
    h3 = c("mu_phi = 0.545 with |phi_i| < 1",rep("",13))
    h4 = c("mu_phi = 0.525 with phi_i in [-1 + epsilon, 1] for some epsilon > 0",rep("",13))
  }
  rv0 = matrix(,1,14)
  h = rbind(h2,h3,h[1:3,],rv0,h[4:6,],rv0,h[7:9,],rv0,h[10:12,],rv0,h4,h[13:15,],rv0,h[16:18,],rv0,h[19:21,],rv0,h[22:24,])
  rownames(h) <- NULL
  colnames(h) <- NULL
  writeData(wb,sn, x = h, startCol = 1, startRow = 3,colNames = FALSE, rowNames = FALSE)
  mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)
  mergeCells(wb, sheet = sn, cols = 4:6, rows = 2)
  mergeCells(wb, sheet = sn, cols = 8:10, rows = 2)
  mergeCells(wb, sheet = sn, cols = 12:14, rows = 2)
  mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 4)
  mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 21)
  
  addStyle(wb,sn,style = right, rows = 5:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = center_style, rows = c(2:3),cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = left, rows = c(4,21),cols = 1:(ncol(h)), gridExpand = T)
  addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:6,8:10,12:14), gridExpand = T,stack=T)
  addStyle(wb,sn,style = bbs, rows = c(3,4,21),cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
  addStyle(wb,sn,style = lbs, rows = (nrow(h)+2),cols = 1:(ncol(h)), gridExpand = T,stack=T)
  rm(tab,tab2)
}

##############################################################################

# Save the workbook to an Excel file
saveWorkbook(wb, file = "HetroAR_MC_FDAC_moments_exp_m1.xlsx",overwrite = TRUE)
cat("The MC results have been written to the excel file HetroAR_MC_FDAC_moments_exp_m1.xlsx.")


## Figure: power functions of FDAC estimator with Gaussian errors without and with GARCH effects
name = "Figure S.4.png"
##########################################################################################################################################################################
#### Define x axis by alternative values
truev = 0.5;
mcm = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01))
intv = 0.2

png(name, units="in", width=28, height=48, res=50)
par(mai=c(2,2,2,2),xpd=TRUE, mfrow=c(4,2),oma=c(6,3,3.5,2)) # (b,l,t,r)

Nlist = c(100,1000,2500,5000); Tlist = c(4,5,6,10)
NT = expand.grid(Nlist,Tlist); NT2 = cbind(NT[,2],NT[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)

#### Create empty matrices to store bias, RMSE, etc.
m1_bias = matrix(,length(NT_l),2)
m1_rmse = matrix(,length(NT_l),2)
m1_size = matrix(,length(NT_l),2)
m1_pw = matrix(,length(NT_l),2*length(mcm))

#### Get the respective simulation results
case_l = c(5,2);
for (j in c(1,2)) {
  case = case_l[j]; init = 1; m0 = 100; GAUSSIAN = 1; GARCH = 0; b = 1; kappa = 2
  fn = paste("exp_fdac_mm_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
  load(fn)
  mb = 0.5; rep = 2000

  for (k in 1:length(NT_l)){
    kid = NT_l[k]
    Tobs = NT[kid,1]; N = NT[kid,2];
    name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
    Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
    s = (Tid-1)*length(Nlist)+Nid;
    if (s<10) {
      name = paste(name,s,sep="_0")
    } else {
      name = paste(name,s,sep="_")
    }
    mc_results1 = get(name);rm(name)
    if (Tobs>3) { # Mean
      results_m1 = mc_results1[1:2000,1]
      results_m1_sd = mc_results1[1:2000,2]
      m1_bias[k,j] = Bias(results_m1, mb)
      m1_rmse[k,j] = RMSE(results_m1, mb)
      e1 = cbind(results_m1,results_m1_sd)
      m1_pw[k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,mcm,rep) # length = no. alternative values
      m1_size[k,j] = m1_pw[k,501*(j-1)+(501-1)/2+1]; rm(e1)
    }
  }
}

power1 = t(m1_pw[,1:501]) # FDAC estimator
power2 = t(m1_pw[,502:1002]) # AAH estimator

#### T=4
d1 = 201; d2= 301
mx = mcm[d1:d2]; xmin = min(mx); xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power1[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis

d1 = 201; d2= 301; mx = mcm[d1:d2]; xmin = min(mx)
xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power2[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis

### T=10
d1 = 201; d2= 301; mx = mcm[d1:d2]; xmin = min(mx)
xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power1[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


d1 = 201; d2= 301; mx = mcm[d1:d2]; xmin = min(mx)
xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power2[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


mtext(expression(paste("Gaussian errors without GARCH effects")), side = 3, cex = 4, font = 2,
      line = -3,
      outer = TRUE)


mtext(expression(paste("Non-Gaussian errors without GARCH effects")), side = 3, cex = 4, font=2,
      line = -180,
      outer = TRUE)


GAUSSIAN = 0;
GARCH = 0;
Nlist = c(100,1000,2500,5000); Tlist2 = c(4,5,6,8,10)
NT = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)

#### Create empty matrices to store bias, RMSE, etc.
m1_bias = matrix(,length(NT_l),2)
m1_rmse = matrix(,length(NT_l),2)
m1_size = matrix(,length(NT_l),2)
m1_pw = matrix(,length(NT_l),2*length(mcm))

case_l = c(5,2);
for (j in c(1,2)) {
  case = case_l[j]; m0 = 100
  fn = paste("exp_fdac_mm_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
  load(fn); rep = 2000
  mb = 0.5

  for (k in 1:length(NT_l)){
    kid = NT_l[k]
    Tobs = NT[kid,1]; N = NT[kid,2];
    name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
    Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
    s = (Tid-1)*length(Nlist)+Nid;
    if (s<10) {
      name = paste(name,s,sep="_0")
    } else {
      name = paste(name,s,sep="_")
    }
    mc_results1 = get(name);rm(name)
    if (Tobs>3) { # Mean
      results_m1 = mc_results1[1:2000,1]
      results_m1_sd = mc_results1[1:2000,2]
      m1_bias[k,j] = Bias(results_m1, mb)
      m1_rmse[k,j] = RMSE(results_m1, mb)
      e1 = cbind(results_m1,results_m1_sd)
      m1_pw[k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,mcm,rep) # length = no. alternative values
      m1_size[k,j] = m1_pw[k,501*(j-1)+(501-1)/2+1]; rm(e1)
    }
  }
}


power1 = t(m1_pw[,1:501])
power2 = t(m1_pw[,502:1002])

#### T=4
d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power1[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power2[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


### T=10
d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power1[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power2[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("n=100", "n=1,000", "n=5,000"), lty = c(4,2,1),col = c(1,2,4), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 5, seg.len=5, bty = 'n')

dev.off()
###################################################################################################################################################################

## Figure: power functions of FDAC estimator with Gaussian and non-Gaussian errors without GARCH effects
name = "Figure S.5.png"
##########################################################################################################################################################################
#### Define x axis by alternative values
truev = 0.5;
mcm = as.matrix(seq(as.numeric(mb-2.5),as.numeric(mb+2.5),by=0.01))
intv = 0.2

png(name, units="in", width=28, height=48, res=50)
par(mai=c(2,2,2,2),xpd=TRUE, mfrow=c(4,2),oma=c(6,3,3.5,2)) # (b,l,t,r)

Nlist = c(100,1000,2500,5000); Tlist2 = c(4,5,6,10)
NT = expand.grid(Nlist,Tlist); NT2 = cbind(NT[,2],NT[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)

#### Create empty matrices to store bias, RMSE, etc.
m1_bias = matrix(,length(NT_l),2)
m1_rmse = matrix(,length(NT_l),2)
m1_size = matrix(,length(NT_l),2)
m1_pw = matrix(,length(NT_l),2*length(mcm))

#### Get the respective simulation results
case_l = c(5,2);
for (j in c(1,2)) {
  case = case_l[j]; init = 1; m0 = 100; GAUSSIAN = 1; GARCH = 0; b = 1; kappa = 2
  fn = paste("exp_fdac_mm_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
  load(fn)
  mb = 0.5; rep = 2000
  
  for (k in 1:length(NT_l)){
    kid = NT_l[k]
    Tobs = NT[kid,1]; N = NT[kid,2];
    name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
    Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
    s = (Tid-1)*length(Nlist)+Nid;
    if (s<10) {
      name = paste(name,s,sep="_0")
    } else {
      name = paste(name,s,sep="_")
    }
    mc_results1 = get(name);rm(name)
    if (Tobs>3) { # Mean
      results_m1 = mc_results1[1:2000,1]
      results_m1_sd = mc_results1[1:2000,2]
      m1_bias[k,j] = Bias(results_m1, mb)
      m1_rmse[k,j] = RMSE(results_m1, mb)
      e1 = cbind(results_m1,results_m1_sd)
      m1_pw[k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,mcm,rep) # length = no. alternative values
      m1_size[k,j] = m1_pw[k,501*(j-1)+(501-1)/2+1]; rm(e1)
    }
  }
}

power1 = t(m1_pw[,1:501]) # FDAC estimator
power2 = t(m1_pw[,502:1002]) # AAH estimator

#### T=4
d1 = 201; d2= 301
mx = mcm[d1:d2]; xmin = min(mx); xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power1[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis

d1 = 201; d2= 301; mx = mcm[d1:d2]; xmin = min(mx)
xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power2[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis

### T=10
d1 = 201; d2= 301; mx = mcm[d1:d2]; xmin = min(mx)
xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power1[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


d1 = 201; d2= 301; mx = mcm[d1:d2]; xmin = min(mx)
xmax = max(mx); ymin = 0; ymax = 1; intv=0.1
plot(mx, power2[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


mtext(expression(paste("Gaussian errors without GARCH effects")), side = 3, cex = 4, font = 2,
      line = -3,
      outer = TRUE)


mtext(expression(paste("Gaussian errors with GARCH effects")), side = 3, cex = 4, font=2,
      line = -180,
      outer = TRUE)


GAUSSIAN = 1;
GARCH = 1;
Nlist = c(100,1000,2500,5000); Tlist2 = c(4,5,6,8,10)
NT = expand.grid(Nlist2,Tlist2); NT2 = cbind(NT2[,2],NT2[,1])
NT_l = c(1,2,4,5,6,8,9,10,12,13,14,16)

#### Create empty matrices to store bias, RMSE, etc.
m1_bias = matrix(,length(NT_l),2)
m1_rmse = matrix(,length(NT_l),2)
m1_size = matrix(,length(NT_l),2)
m1_pw = matrix(,length(NT_l),2*length(mcm))

case_l = c(5,2);
for (j in c(1,2)) {
  case = case_l[j]; m0 = 100
  fn = paste("exp_fdac_mm_c",case,"_i",init,"_m",m0,"_b",b,"_kp",kappa,"_",GAUSSIAN,GARCH,".RData",sep="")
  load(fn); rep = 2000
  mb = 0.5
  
  for (k in 1:length(NT_l)){
    kid = NT_l[k]
    Tobs = NT[kid,1]; N = NT[kid,2];
    name = paste("mm","_c",case,"_m",m0,"_b",b,sep="")
    Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
    s = (Tid-1)*length(Nlist)+Nid;
    if (s<10) {
      name = paste(name,s,sep="_0")
    } else {
      name = paste(name,s,sep="_")
    }
    mc_results1 = get(name);rm(name)
    if (Tobs>3) { # Mean
      results_m1 = mc_results1[1:2000,1]
      results_m1_sd = mc_results1[1:2000,2]
      m1_bias[k,j] = Bias(results_m1, mb)
      m1_rmse[k,j] = RMSE(results_m1, mb)
      e1 = cbind(results_m1,results_m1_sd)
      m1_pw[k,(501*(j-1)+1):(501*(j-1)+501)] = Power(mb,e1,mcm,rep) # length = no. alternative values
      m1_size[k,j] = m1_pw[k,501*(j-1)+(501-1)/2+1]; rm(e1)
    }
  }
}


power1 = t(m1_pw[,1:501])
power2 = t(m1_pw[,502:1002])

#### T=4
d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power1[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power2[d1:d2,1], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=4")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,2], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,3], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


### T=10
d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power1[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Homogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power1[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power1[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


d1 = 201; d2= 301
mx = mcm[d1:d2]
xmin = min(mx)
xmax = max(mx)
ymin = 0
ymax = 1
intv=0.1
plot(mx, power2[d1:d2,10], type="l", col=1, lwd=5, lty=4,xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs="i",yaxs="i",ylab="",xlab="", main=expression(paste("Heterogeneous with T=10")),cex.lab=4,  axes = F, cex.main=4.5, cex.sub=4)
lines(mx, power2[d1:d2,11], type="l", col=2, lwd=5, lty=2,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
lines(mx, power2[d1:d2,12], type="l", col=4, lwd=5, lty=1,xlim=c(xmin, xmax), ylim=c(ymin, ymax))
axis(side=1, at=round(seq(xmin, xmax, by=intv),1), mgp=c(2, 3.5, 0), cex.axis=3.5)
axis(side=2, at=seq(0, 1, by=intv), mgp=c(4, 2, 0), cex.axis=3.5)######## Visualization
arrows(truev,-0.01, truev, 1, length = 0,lty =3,lwd=3) # denote the true value
arrows(xmin,0.05,xmax,0.05, length = 0,lty = 3, lwd=3) # denote a line of 5% probability
mtext("power", side = 2, line = 5, at = 1 , cex = 3)
mtext("5%", side = 1, line = -5, at = xmax , cex = 2.5)
mtext(expression(mu[phi]), side = 1, line = 8, at = xmax , cex = 3) # label of x-axis
mtext(expression(mu[paste(phi,",0")]), side = 1, line = 8, at = truev , cex = 3) # label of x-axis


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("n=100", "n=1,000", "n=5,000"), lty = c(4,2,1),col = c(1,2,4), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 5, seg.len=5, bty = 'n')

dev.off()
###################################################################################################################################################################





