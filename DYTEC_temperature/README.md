[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **DYTEC_temperature** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿Name of Quantlet: 'DYTEC_temperature'

Published in: Dynamcis of Tail Event Curves

Description: 'This is the application of  DYTEC algorithms to the temperature data'

Keywords: expectile, dynamics, factor model


Author: Petra Burdejova

Submitted:  Nov 28 2017 by Petra Burdejova

Input: 'temperature.txt'

Output:  Compute bata-matrix and plot  compare forecast for expectile curve
```

### R Code
```r

# ---------------------------------------------------------------------
# Paper:       	P. Burdejova and W.K. Härdle
#			"Dynamics of Tail Event Curves"
# ---------------------------------------------------------------------
# Quantlet:    	DYTEC_temperature
# ---------------------------------------------------------------------
# Description: 	Do the simulation for DYTEC algorithm
# ---------------------------------------------------------------------
# Author:      	Petra Burdejova
# ---------------------------------------------------------------------

library(expectreg)
library(fda)
library(matrixcalc)
library(gglasso)
library(doParallel)
#library(parallel)
library(foreach)
library(fBasics)

install.packages("vars")
library("vars")


wgt_fun <- function(tau,basis,mtx_hat,y_expl){	
  wgt_vec <- y_expl < (basis %*% mtx_hat)	
  wgt_vec <- abs(wgt_vec - tau)
  return(wgt_vec)
}

my_tau=0.5
tau=0.5

#stat=1

for (stat in 1:1){
  data0 = t(temp_list_stations[[stat]])
  data0= data0 #[1:20,30:150]
  data0_exp = t(exp_list_stations[[stat]])
  data0_exp= data0_exp #[1:20,30:150]
  fcast_dytec = dytec_fcast(data0) 
  fcast_var = var_fcast(data0_exp)
} # end of for cycle stat



# 
# plot(data0[53,])
# lines(data0_exp[53,], col="red")
# lines(fcast_dytec)
# lines(fcast_var)
#lines(fcast_dytec_sm)
lines(unlist(fcast_dytec_sm$values), col="blue")

cl<-makeCluster(2)
registerDoParallel(cl)

sim_output <- 
  foreach (stat=c(1:159), .packages=c('fda','gglasso','expectreg','matrixcalc','vars')) %dopar% {
    data0 = t(temp_list_stations[[stat]])
    data0= data0 #[1:20,30:150]
    data0_exp = t(exp_list_stations[[stat]])
    data0_exp= data0_exp #[1:20,30:150]
    fcast_dytec = dytec_fcast(data0) 
    fcast_dytec_sm <- expectreg.ls(fcast_dytec~rb(c(1:365),type="pspline"), smooth="aic",expectile=0.5)
    fcast_dytec_sm_disc <- unlist(fcast_dytec_sm$values)
    fcast_var = var_fcast(data0_exp)
    compare_fcast(fcast_dytec_sm_disc,fcast_var,data0_exp[53,])
  }
stopCluster(cl)


compare_fcast <- function(dytec_fun, var_fun, exp_fun){
  p=length(exp_fun)
  mse_dytec <- sum(abs(exp_fun-dytec_fun)^2)/p  
  mse_var <- sum(abs(exp_fun-var_fun)^2)/p
  return(list(mse_dytec,mse_var) )
}

# mse_var <- sum(abs(data0_exp[53,]-fcast_var)^2)/365
# mse_dytec <- sum(abs(data0_exp[53,]-fcast_dytec)^2)/365
# mse_var
# mse_dytec


var_fcast <- function(data_exp){
  d1=dim(data_exp)[1]
  d2=dim(data_exp)[2] 
  pca <- prcomp(data_exp[1:(d1-1),])
  components<- pca$rotation[,1:3]
  #dim(components)
  projection <- pca$x[,1:3]
  #dim(projection)
  mean <- apply(data_exp[1:(d1-1),], 2,mean)

    # var na y tj. projection
    scores_var <- VAR(projection, p = 3)
    scores_f <-(predict(scores_var,n.ahead=1, ci=0))
    # forecast(a_pred, h=1)$mean }
    LOC<-c("PC1","PC2","PC3")
    scores_fcast <- sapply(scores_f $fcst[LOC], function (k) k[ , 1])
    e_fcast_var <- components%*% scores_fcast
    return(e_fcast_var+mean)
} # end of var_fcast function




dytec_fcast <- function(data){
  d1=dim(data)[1]
  d2=dim(data)[2] 
  nbas_time=round(d1+1)
  nbas_space=round(d2/2)
  
  #--- create time basis
  time_basis <- dytec_create_time_basis(d1,nbas_time)
  time_basis_matrix <- time_basis[[1]]
  # if pca used, then nbas_space could change
  nbas_time <- time_basis[[2]]
  
  space_basis <- dytec_create_space_basis(data,opt=1,nbas_space)
  space_basis_matrix <- space_basis[[1]]
  # if pca used, then nbas_space could change
  nbas_space <- space_basis[[2]]
  space_basis_matrix <- t(space_basis_matrix)
  
  # delete last year for modelling
  x_vec <- kronecker(t(space_basis_matrix), time_basis_matrix[1:d1-1,]) # delete last year for modelling
  #dim(x_vec)   # 200 x 1400
  data_vec <- vec(data[1:d1-1,])   # delete last year for modelling
  #dim(data_vec)
  
  #--- minimize for mean
  m1 <- gglasso(x=x_vec, y=data_vec, intercept=FALSE)
  #dim(m1$beta)
  beta_hat_mtx <- matrix(m1$beta[,100],nrow=nbas_time, ncol=nbas_space)
  
  J=dim(x_vec)[2]
  #my_tau=0.8
  wgt_old<- wgt_fun(my_tau,x_vec,m1$beta[,100], data_vec)
  w_iter_k <- 0
  min_change <- length(wgt_old)
  repeat {
    w_iter_k= w_iter_k+1
    print(paste("tau=", my_tau ,"iter:",w_iter_k))
    wgt_x_vec <- x_vec
    for ( j in 1:J)  {wgt_x_vec[,j] <- wgt_old * x_vec[,j] }
    mk <- gglasso(x= wgt_x_vec , y= wgt_old * data_vec, intercept=FALSE)
    beta_k <- mk$beta[,100]
    wgt_new <- wgt_fun(my_tau, x_vec, beta_k, data_vec )
    print(paste("diff.w.:",sum(wgt_old != wgt_new)))
    # avoid cycling
    if (sum(wgt_old != wgt_new)>=min_change){break}
    if (sum(wgt_old != wgt_new)<min_change) {min_change <- sum(wgt_old != wgt_new)}
    wgt_old <- wgt_new
  } # end of repeat
  beta_hat_mtx_k <- matrix(beta_k,nrow=nbas_time, ncol=nbas_space)
  y_fit <-  (time_basis_matrix%*% beta_hat_mtx_k %*% space_basis_matrix)
  return(y_fit[d1,]) #return just forecast (made with d1-1 curves)
} # end of dytec fcast

  


  
  

```

automatically created on 2018-09-04