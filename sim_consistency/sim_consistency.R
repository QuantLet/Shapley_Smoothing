library(parallel)


ISE_fct = function(m, obs){

  library(np)
  library(pracma)
  library(cubature)
  library(simstudy)
  library(MASS)
  
  source("functions.R") 
  
  cova<<-0
  sigma_sim<<-matrix(c(4, cova, cova,
                   cova, 4, cova,
                   cova, cova, 4), nrow=3, ncol=3)

  source("integral_population.R")
  source("integral_estimation.R")
  source("shapley_int.R")
  source("SE_vec_int.R")
  
  
  g1 <<- function(X){ return( -sin(2*X[,1]) ) } 
  g2 <<- function(X){ return( cos(3*X[,2])   ) } 
  g3 <<- function(X){ return( 0.5*X[,3] ) } 
  int = function(X){
    x1 = X[,1]
    x2 = X[,2]
    return( 2*cos(x1)*sin(2*x2) ) 
  }
  
  
  l <<- -2; u <<- 2; N<<- obs
  l_int <<- l; u_int <<- u
  d <<- 3
  
  
  X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))
  
  #DGP
  Y <<- g1(X) + g2(X) + g3(X) + int(X) + rt(n=nrow(X), df=5)
  
  #All possible subsets
  subs <<- subsets(X)
  
  #Get model fits and sort them in a list
  model_list <<- model_list_fct(subs=subs) 
  
  
  # Component-based
  ISE_res1=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)
  ISE_res2=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)
  ISE_res3=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)
  
  ISE1 = ISE_res1$integral 
  ISE2 = ISE_res2$integral 
  ISE3 = ISE_res3$integral 
  

  
  
  
  # Integral-based
  ISE_res1_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)
  ISE_res2_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)
  ISE_res3_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)
  
  ISE1_int = ISE_res1_int$integral 
  ISE2_int = ISE_res2_int$integral  
  ISE3_int = ISE_res3_int$integral
  
  
  return(c(ISE1, ISE1_int, ISE2, ISE2_int, ISE3, ISE3_int))
}

res1=mclapply(1:6000, ISE_fct, mc.cores=40, obs=300)
res2=mclapply(1:6000, ISE_fct, mc.cores=40, obs=500)
res3=mclapply(1:6000, ISE_fct, mc.cores=40, obs=1000)
res4=mclapply(1:6000, ISE_fct, mc.cores=40, obs=2000)


results1=matrix(unlist(res1), byrow = FALSE, ncol=6000)
results2=matrix(unlist(res2), byrow = FALSE, ncol=6000)
results3=matrix(unlist(res3), byrow = FALSE, ncol=6000)
results4=matrix(unlist(res4), byrow = FALSE, ncol=6000)

a=rowMeans(results1)
b=rowMeans(results2)
c=rowMeans(results3)
d=rowMeans(results4)

tab_final=rbind(a,b,c,d)
