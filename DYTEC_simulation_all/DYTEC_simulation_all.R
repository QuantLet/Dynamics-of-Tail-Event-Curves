# ---------------------------------------------------------------------
# Paper:       	P. Burdejova and W.K. Härdle
#			"Dynamics of Tail Event Curves"
# ---------------------------------------------------------------------
# Quantlet:    	DYTEC_simulation_all
# ---------------------------------------------------------------------
# Description: 	Do the simulation for DYTEC algorithm
# ---------------------------------------------------------------------
# Author:      	Petra Burdejova
# ---------------------------------------------------------------------


#------------------parameters
getwd()
setwd("/Users/PetraB/GitHub/DSFM_for_expectiles/00_DYTEC_CODE/00_NEW")
params<-read.table("params.txt")# file contain parameters of the simulation (sample size, taus)

K=200 #number of simulation turns
#run with 2 components
nc = 2

#optim params
tol = 1e-10
max.iter=30


#----------------mean and component functions
m<-function(t){(mu=1+t+exp(-(t-0.6)^2/0.05))
  return(mu)
}
ff<-function(t){f1=sqrt(2)*sin(2*pi*t)
f2=sqrt(2)*cos(2*pi*t)	
return(cbind(f1,f2))
}


#--- FUNCTION: simulates curves 
simcurv <- function
###X_ij =mu(t_j)+f_1(t_j)*a_i1+f_2(t_j)*a_i2+...+e_ij, a_i1~N(0,s_1), a_i2~N(0,s_2)
###with option (1) e_ij~ N(0,se), (2) e_ij~t(se), (3)e_ij~N(0,mu(t_j)^2*se), 
### (4)  e_ij~ logN(0,se), (5)  e_ij~ U(0,se)+U(0,se/2)
###Returns vectorized X and X_e theoretical expectile of X stored in a (2pn+4)xK*length(tauseq).
###Rows: vectorized X and X_e resp.+ 4 parameters: tau level,scenario number,n,p. Columns: K simulation runs for each tau in tauseq.
###The default option returns X with the normal error.
(mu, 
 ###mean curve, px1, p is the data dimension
 f, 
 ###pxk for k principal directions
 option = 1, 
 ###option for the error distribution
 s,
 #c(s1, s2,.., se) parameters of the distributions for a and e
 tauseq,
 ###the expectile level, 
 ###s1,s2,se are variances for normal distributions, se are degrees of freedom for t distr.
 N,P,K
 ### dimensions of the data: N number of curves, P grid, K number of simulation runs
){
  set.seed(1234567)
  k=length(s)-1
  se=s[k+1]
  mua=array(mu, c(P,N,K)); fa=array(NA, c(P,N,K))
  a=sweep(array(rnorm(k*N*K),c(k,N,K)),MARGIN=1,sqrt(s[1:k]),`*`)
  output=matrix(NA,(2*N*P+4),length(tauseq)*K)
  for (y in 1:length(tauseq)){
    set.seed(123)
    tau=tauseq[y]
    X  = array(NA,c(P,N,K))
    X_e= array(NA,c(P,N,K)) #the individual expectile
    if (option==1){
      r=array(rnorm(N*P*K),c(P,N,K))
      e=sqrt(se)*r
      etau=array(enorm(tau, sd=sqrt(se)), c(P,N,K))
    } else if (option==2) {
      e=array(rt(N*P*K,se),c(P,N,K))
      etau=array(et(tau, se), c(P,N,K))
    } else if (option==3) {
      r=array(rnorm(N*P*K),c(P,N,K))
      e=sweep(r,MARGIN=1,sqrt(se*mu),`*`)
      etau=numeric(P)
      for (i in 1:P){
        etau[i]=enorm(tau, sd=sqrt(se*mu[i]))
      }
      etau=array(etau,c(P,N,K))
    } else if (option==4) {
      e=rlnorm(N*P*K,sdlog=sqrt(se)); e=array(e,c(P,N,K))
      etau=array(elnorm(tau, sdlog=sqrt(se)), c(P,N,K))
    } else if (option==5) {
      ee=runif(100000,max=se)+runif(100000,max=se)
      e=runif(N*P*K,max=se)+runif(N*P*K,max=se); e=array(e,c(P,N,K))
      etau=array(expectile(matrix(ee,1,100000),alpha=(tau-0.5)),c(P,N,K))
    } else {stop("'option' takes integer values 1 to 5")}	
    for (j in 1:K){
      fa[,,j]=f%*%a[,,j]
    }
    X  = mua+fa+e
    X_e= mua+fa+etau
    output[,(K*(y-1)+1):(K*y)]=rbind(matrix(X,N*P,K),matrix(X_e,N*P,K),rep(tau,times=K),rep(option,times=K),rep(n,times=K),rep(p,times=K))
    #output[,(K*(y-1)+1):(K*y)]=rbind(matrix(X,N*P,K),matrix(X_e,N*P,K))
  }
  return(output)
}


#params0= params
#params = params0[1,]
#params[3]=6
#params[4]=1

noptions=4 #number of scenarios

my_tau=0.8


#----------------generate the curves (takes some time: 15 min)
XX=as.list(1:(nrow(params)))
for (q in 1:nrow(params)){
  cat("q",q, "\n")
  n=params[q,1]; p=params[q,2]
  #variance parameter of N(0,s)	
  t=seq(0,1,length=p)
  s1=params[q,3]
  s2=params[q,4] 
  se=params[q,5]
  s=c(s1,s2,se)
  #degrees of freedom for t
  df=params[q,6]
  mu=m(t); f=ff(t)
  #X=NULL
  X=as.list(1:noptions)
  for (h in 1:noptions){
    cat("q",q, "h", h, "\n")
    option=h; s[3]=ifelse(option==2,df,se)
    #X=cbind(X,simcurv(mu=mu,f=f,option=h,s=s,tauseq=my_tau, N=n,P=p,K=K))#output matrix 2np x length(tauseq)K
    X[[h]]=simcurv(mu=mu,f=f,option=h,s=s,tauseq=my_tau, N=n,P=p,K=K)#output matrix 2np x length(tauseq)K
  }
  XX[[q]]<-X
}

### XX
###vectorized matrix of dimensions (2pn+4),1) x length(tauseq)K.
##first pn elements contain simulated curves (data),
##next pn elements contain theoretical expectile curves, 
##last elements contains tau expectile level,scenario, n,p.

#scenario=1
#eps_option=1

avg_MSE_all <- matrix(nrow=6, ncol=4)


for (scenario in c(1:6)){
  for (eps_option in c(1:4)){
    X_sim= XX[[scenario]][[eps_option]]
    p=X_sim[dim(X_sim)[1],1]
    n=X_sim[dim(X_sim)[1]-1,1]
    # eps_option X_sim[dim(X_sim)[1]-2,1]
    tau=X_sim[dim(X_sim)[1]-3,1]
    #ind_sim=1 
    MSE=rep(0,K)
    for (ind_sim in c(1:K)){
      Xd  = matrix(X_sim[1:(n*p),ind_sim],,n) #extract data
      X_e  = matrix(X_sim[(n*p+1):(2*n*p),ind_sim],p,n) #extract theoretical expectiles
      X_fit <- dytec_estim(Xd,tau)
      MSE[ind_sim] = sum((t(X_e) - X_fit)^2)/(n*p)
    }
    avg_MSE_all[scenario][eps_option] <- sum(MSE)/K
    
    
  } # end of eps_option cycle
} # end of eps_option cycle

data <-Xd

matplot(Xd,type="l")

wgt_fun <- function(tau,basis,mtx_hat,y_expl){	
  wgt_vec <- y_expl < (basis %*% mtx_hat)	
  wgt_vec <- abs(wgt_vec - tau)
  return(wgt_vec)
}


dytec_estim <- function(data,my_tau){
  
  data <- t(data)
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
  
  #dim(space_basis_matrix)
  #dim(time_basis_matrix)
  
  x_vec <- kronecker(t(space_basis_matrix), time_basis_matrix) # delete last year for modelling
  #dim(x_vec)   # 200 x 1400
  data_vec <- vec(data) 
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
  return(y_fit)
  
} # end of function est_dytec




