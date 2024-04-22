###### Pesaran and Yang (2024) "Hetro. AR panels"
rm(list=ls())

## Please change the path to the directory where the simulation outcomes are stores. 
# setwd('~/Downloads') 

## Please install the following packages. 
# list.of.packages <- c("openxlsx")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(openxlsx)


## Compare FDAC and MSW estimators
##############################################################################
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

pwlb=-2.5; pwrb=2.5; pwintv=0.01;
pwgrid = length(as.matrix(seq(as.numeric(0+pwlb),as.numeric(0+pwrb),by=pwintv)))
Power = function(truev,estsd,rep){
  h1 = as.matrix(seq(as.numeric(truev+pwlb),as.numeric(truev+pwrb),by=pwintv)) # values of alternatives (v0-0.5,v0+0.5) length = 101
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

# Create a style for center-aligning the numbers
##############################################################################
center_style <- createStyle(halign = "center")
right <- createStyle(halign = "right")
left <- createStyle(halign = "left")
bbs <- createStyle(border = "bottom",borderStyle="thin")
lbs <- createStyle(border = "bottom",borderStyle="double")
tbs <- createStyle(border = "top",borderStyle="double")
wrap_text_style <- createStyle(wrapText = TRUE)

pwlb=-2.5; pwrb=2.5; pwintv=0.01;
pwgrid = length(as.matrix(seq(as.numeric(0+pwlb),as.numeric(0+pwrb),by=pwintv)))
k=2;
##############################################################################

wb <- createWorkbook()

## Table S.11
##############################################################################
# Load simulation results
case_l = c(1,2); init = 1; m0 = 100; b = 1; kappa = 2; mb_l = c(0.4, 0.5)
fn1 = paste("exp_fdac_msw_exp_c",case_l[1],"_i",init,"_m",m0,"_b",b,"_kp",kappa,".RData",sep="")
fn2 = paste("exp_fdac_msw_exp_c",case_l[2],"_i",init,"_m",m0,"_b",b,"_kp",kappa,".RData",sep="")
en = 1; temp1 = tryCatch(load(fn1),error=function(e) print(en))
temp2 = tryCatch(load(fn2),error=function(e) print(en))

if ( temp1[1]!=1 & temp2[1]!=1) {
  # Create a new workbook
  sn = "Table S.11"
  addWorksheet(wb, sn)
  writeData(wb, sn, x = "Table S.11: Bias, RMSE, and size of FDAC and MSW estimators of mu_phi = E(phi_i) in heterogeneous panel AR(1) models with uniformly distributed phi_i and Gaussian errors without GARCH effects", startCol = 1, startRow = 1, colNames = FALSE)
  writeData(wb, sn, x = "mu_phi = 0.4 with |phi_i|<1", startCol = 4, startRow = 2, colNames = FALSE)
  writeData(wb, sn, x = "mu_phi = 0.5 with phi_i in [-1 + epsilon, 1] for some epsilon > 0", startCol = 13, startRow = 2, colNames = FALSE)
  
  Nlist = c(100,1000); Tlist = c(4,6,10)
  NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
  
  m1_bias = matrix(,2*dim(NT)[1],2)
  m1_rmse = matrix(,2*dim(NT)[1],2) 
  m1_size = matrix(,2*dim(NT)[1],2)
  m1_pw = matrix(,2*dim(NT)[1],2*pwgrid)
  
  for (case in case_l) {
    mb = mb_l[case]
    #### Get the respective simulation results
    for (k in 1:dim(NT)[1]){ 
      Tobs = NT[k,1]; N = NT[k,2]
      name = paste("mc","_c",case,"_m",m0,sep="")
      Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
      s = (Tid-1)*length(Nlist)+Nid;
      if (s<10) {
        name = paste(name,s,sep="_0")
      } else {
        name = paste(name,s,sep="_") 
      }
      mc_results1 = get(name);rm(name)
      
      for (j in 1:2) { # 2 estimators: FDAC and MSW
        #### Calculate bias, RMSE, size and powers for moments
        idl = which(is.na(mc_results1[(2-j)*4+1,])==0); rep = length(idl)
        results_m1 = mc_results1[(2-j)*4+1,idl] # row = 1, 5
        results_m1_sd = mc_results1[(2-j)*4+2,idl] # row = 2, 6
        m1_bias[(case-1)*dim(NT)[1]+k,j] = Bias(results_m1, mb)
        m1_rmse[(case-1)*dim(NT)[1]+k,j] = RMSE(results_m1, mb)
        e1 = cbind(results_m1,results_m1_sd)
        m1_pw[(case-1)*dim(NT)[1]+k,(pwgrid*(j-1)+1):(pwgrid*(j-1)+pwgrid)] = Power(mb,e1,rep) # length = no. alternative values
        m1_size[(case-1)*dim(NT)[1]+k,j] = m1_pw[(case-1)*dim(NT)[1]+k,pwgrid*(j-1)+(pwgrid-1)/2+1]; rm(e1) 
      }
    }
  }
  
  #### Generate Table with empty entries
  #### Set Table Size
  #########################################################################
  nr2 = dim(NT)[1]; 
  Tlist2 = c("4","6","10"); Nlist2 = c('100','1,000','100','1,000',"100","1,000")
  nc = 2
  tab = matrix(,nr2,2+(nc+1)*3*2-1) 
  tab2 = matrix(,nr2,2+(nc+1)*3*2-1) 
  #########################################################################
  
  #### Put numbers into the respective cells
  #########################################################################
  for (r2 in 1:nr2) {
    tab2[r2,1] = NT[r2,1]
    tab2[r2,2] = Nlist2[r2]    
    tab[r2,(2+1):(2+nc)] = m1_bias[r2,]
    tab[r2,(2+nc+1+1):(2+nc+1+nc)] = m1_rmse[r2,]
    tab[r2,(2+nc+1+nc+1+1):(2+nc+1+nc+1+nc)] = m1_size[r2,]
    
    tab[r2,(5+3*nc+1):(5+3*nc+nc)] = m1_bias[r2+nr2*1,]
    tab[r2,(5+3*nc+nc+1+1):(5+3*nc+nc+1+nc)] = m1_rmse[r2+nr2*1,]
    tab[r2,(5+3*nc+nc+1+nc+1+1):(5+3*nc+nc+1+nc+1+nc)] = m1_size[r2+nr2*1,]
  }
  #########################################################################
  
  #### Change the formats of each column
  #########################################################################
  for (col in c(3:4,6:7,12:13,15:16)) {   # bias, RMSE
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
  for (col in c(9:10,18:19)) {   # size (*100)
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
  
  cv0 = matrix(,nrow = nrow(tab2),1)
  h = cbind(tab2[,1:2],cv0,tab2[,3:ncol(tab2)])
  h1 = c("","","","Bias","","","RMSE","","","Size (*100)","","","Bias","","","RMSE","","","Size (*100)","")
  h2 = c("T","n","","FDAC","MSW","","FDAC","MSW","","FDAC","MSW","","FDAC","MSW","","FDAC","MSW","","FDAC","MSW")
  rv0 = matrix(,1,ncol(h))
  h = rbind(h1,h2,h[1:2,],rv0,h[3:4,],rv0,h[5:6,])
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
  
  addStyle(wb,sn,style = center_style, rows = 2:4,cols = 1:(ncol(h)), gridExpand = TRUE)
  addStyle(wb,sn,style = right, rows = 5:(nrow(h)+2),cols = 1:(ncol(h)), gridExpand = TRUE)
  addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = TRUE,stack=TRUE)
  addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:11,13:20), gridExpand = TRUE,stack=TRUE)
  addStyle(wb,sn,style = bbs, rows = 3,cols = c(4:5,7:8,10:11,13:14,16:17,19:20), gridExpand = TRUE,stack=TRUE)
  addStyle(wb,sn,style = bbs, rows = 4,cols = 1:(ncol(h)), gridExpand = TRUE,stack=T)  
  addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = TRUE,stack=TRUE)
  addStyle(wb,sn,style = lbs, rows = (nrow(h)+2),cols = 1:(ncol(h)), gridExpand = TRUE,stack=TRUE)
} 
##############################################################################

## Table S.17
##############################################################################
# Load simulation results
case_list = c(1,2,5); m0_list = c(100,1)
init = 1; b = 1; kappa = 2; mb_l = c(0.4, 0.5, 0.5)
ncase = length(case_list); nm = length(m0_list)
nest = 2

sn = "Table S.17"
addWorksheet(wb, sn)
writeData(wb, sn, x = "Table S.17: Bias, RMSE, and size of FDAC and MSW estimators of mu_phi = E(phi_i) in heterogeneous and homogeneous panel AR(1) models with Gaussian errors without GARCH effects and different initializations", startCol = 1, startRow = 1, colNames = FALSE)

Nlist = c(100,1000); Tlist = c(4,6,10)
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
m1_bias = matrix(,ncase*dim(NT)[1],nest*nm)
m1_rmse = matrix(,ncase*dim(NT)[1],nest*nm) 
m1_size = matrix(,ncase*dim(NT)[1],nest*nm)
m1_pw = matrix(,ncase*dim(NT)[1],nest*nm*pwgrid)

for (cid in 1:ncase) {
  case = case_list[cid]
  fn1 = paste("exp_fdac_msw_exp_c",case,"_i",init,"_m",m0_list[1],"_b",b,"_kp",kappa,".RData",sep="")
  fn2 = paste("exp_fdac_msw_exp_c",case,"_i",init,"_m",m0_list[2],"_b",b,"_kp",kappa,".RData",sep="")
  en = 1; temp1 = tryCatch(load(fn1),error=function(e) print(en))
  temp2 = tryCatch(load(fn2),error=function(e) print(en))
  mb_l = c(0.4, 0.5, 0.5)
  mb = mb_l[cid]
  
  if ( temp1[1]!=1 & temp2[1]!=1) {
    Nlist = c(100,1000); Tlist = c(4,6,10)
    NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
    
    #### Get the respective simulation results
    for (mid in 1:nm) {
      m0 = m0_list[mid]
      for (k in 1:dim(NT)[1]){ 
        Tobs = NT[k,1]; N = NT[k,2]
        name = paste("mc","_c",case,"_m",m0,sep="")
        Nid = which(Nlist==N); Tid = which(Tlist==Tobs)
        s = (Tid-1)*length(Nlist)+Nid;
        if (s<10) {
          name = paste(name,s,sep="_0")
        } else {
          name = paste(name,s,sep="_") 
        }
        mc_results1 = get(name);rm(name)
        
        for (j in 1:2) { # 2 estimators: FDAC and MSW
          #### Calculate bias, RMSE, size and powers for moments
          idl = which(is.na(mc_results1[(2-j)*4+1,])==0); rep = length(idl)
          results_m1 = mc_results1[(2-j)*4+1,idl] # row = 1, 5
          results_m1_sd = mc_results1[(2-j)*4+2,idl] # row = 2, 6
          
          m1_bias[(cid-1)*dim(NT)[1]+k,mid+2*(j-1)] = Bias(results_m1, mb)
          m1_rmse[(cid-1)*dim(NT)[1]+k,mid+2*(j-1)] = RMSE(results_m1, mb)
          e1 = cbind(results_m1,results_m1_sd)
          m1_pw[(cid-1)*dim(NT)[1]+k,(pwgrid*(mid+2*(j-1)-1)+1):(pwgrid*(mid+2*(j-1)-1)+pwgrid)] = Power(mb,e1,rep) # length = no. alternative values
          m1_size[(cid-1)*dim(NT)[1]+k,mid+2*(j-1)] = m1_pw[(cid-1)*dim(NT)[1]+k,pwgrid*(mid+2*(j-1)-1)+(pwgrid-1)/2+1]; rm(e1)
        }
      }
    }
  }
}

#### Generate Table with empty entries
#### Set Table Size
#########################################################################
Nlist = c(100,1000); Tlist = c(4,6,10); ncase = 3
NT = expand.grid(Nlist,Tlist); NT = cbind(NT[,2],NT[,1])
nc = (nm+1)*nest*3 # bias, RMSE, size
nr = length(Tlist)*(length(Nlist)+1)+1
Tlist2 = c("4","6","10"); Nlist2 = c('100','1,000','100','1,000','100','1,000')
tab = matrix(,ncase*nr,nc+2) 
tab2 = matrix(,ncase*nr,nc+2) 
#########################################################################

#### Put numbers into the respective cells
#########################################################################
N2 = c('100','1,000'); T2 = c(4,6,10)
NT2 = rbind(expand.grid(N2,T2),expand.grid(N2,T2),expand.grid(N2,T2))
cv0 = matrix(,nrow(m1_bias),1)
tab0 = cbind(NT2,cv0,m1_bias[,1:2],cv0,m1_bias[,3:4],cv0,m1_rmse[,1:2],cv0,m1_rmse[,3:4],cv0,m1_size[,1:2],cv0,m1_size[,3:4])
rv0 = matrix(,1,ncol(tab0))
colnames(rv0) = colnames(tab0)
tab = tab0[1:2,]
for (r1 in 2:(3*length(Tlist))) {
  tab = rbind(tab,rv0,tab0[(1+(r1-1)*2):(r1*2),])
}
tab1 = rbind(rv0,tab[1:9,],rv0,tab[10:18,],rv0,tab[19:26,])
tab2 = tab1
#########################################################################

#### Change the formats of each column
#########################################################################
for (col in c(4:5,7:8,10:11,13:14)) {   # bias, RMSE
  v0 = tab1[,col]; v1 = v0
  for (i in 1:length(v0)) {
    if (is.na(v0[i]) == 0) {
      v1[i] = format(round(v0[i], 3), nsmall = 3) 
    } else {
      v1[i] = ""
    }
  }
  tab2[,col] = v1; rm(v0,v1)
}
for (col in c(16:17, 19:20)) {   # size (*100)
  v0 = tab1[,col]; v1 = v0
  for (i in 1:length(v0)) {
    if (is.na(v0[i]) == 0) {
      v1[i] = format(round(v0[i]*100, 1), nsmall = 1) 
    } else {
      v1[i] = ""
    }
  }
  tab2[,col] = v1; rm(v0,v1)
}

h1 = c("","","","Bias",rep("",4),"","RMSE",rep("",4),"","Size (*100)",rep("",4))
h2 = c("","",rep(c("","FDAC","","","MSW",""),3))
h3 = c("T","n/M_0", rep(c("","100","1"),6))
tab2[1,1] = "mu_phi = 0.4 with uniformly distributed |phi_i|<1 for all i"
tab2[11,1] = "mu_phi = 0.5 with uniformly distributed phi_i in [-1 + epsilon, 1] for some epsilon > 0 and all i"
tab2[21,1] = "mu_phi = 0.5 for all i"
h = rbind(h1,h2,h3,tab2)
rownames(h) <- NULL
colnames(h) <- NULL
writeData(wb,sn, x = h, startCol = 1, startRow = 2,colNames = FALSE, rowNames = FALSE)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 1)  
mergeCells(wb, sheet = sn, cols = 4:8, rows = 2)
mergeCells(wb, sheet = sn, cols = 10:14, rows = 2)
mergeCells(wb, sheet = sn, cols = 16:20, rows = 2)
mergeCells(wb, sheet = sn, cols = 4:5, rows = 3)
mergeCells(wb, sheet = sn, cols = 7:8, rows = 3)
mergeCells(wb, sheet = sn, cols = 10:11, rows = 3)
mergeCells(wb, sheet = sn, cols = 13:14, rows = 3)
mergeCells(wb, sheet = sn, cols = 16:17, rows = 3)
mergeCells(wb, sheet = sn, cols = 19:20, rows = 3)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 5)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 15)
mergeCells(wb, sheet = sn, cols = 1:ncol(h), rows = 25)

addStyle(wb,sn,style = center_style, rows = 2:4,cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = left, rows = c(5,15,25),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = right, rows = c(6:14,16:24,26:34),cols = 1:(ncol(h)), gridExpand = T)
addStyle(wb,sn,style = wrap_text_style, rows = 1,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 2,cols = c(4:8, 10:14, 16:20), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = 3,cols = c(4:5,7:8,10:11,13:14,16:17,19:20), gridExpand = T,stack=T)
addStyle(wb,sn,style = bbs, rows = c(4,5,15,25),cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = tbs, rows = 2,cols = 1:(ncol(h)), gridExpand = T,stack=T)
addStyle(wb,sn,style = lbs, rows = (nrow(h)+1),cols = 1:(ncol(h)), gridExpand = T,stack=T)
##############################################################################

# Save the workbook to an Excel file
saveWorkbook(wb, file = "HetroAR_MC_FDAC_vs_MSW_results.xlsx",overwrite = TRUE)
cat("The MC results have been written to the excel file HetroAR_MC_FDAC_vs_MSW_results.xlsx.")



