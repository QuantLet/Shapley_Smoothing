
bootstrap = function(j, N=2000, B=10000){

library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)

source("functions.R") 

cov=0
sigma_sim<<-matrix(c(4, cov,
               cov, 4), nrow=2, ncol=2)


source("integral_population_2d_int.R")
source("integral_estimation_2d.R") 

source("shapley_int.R")

g1 = function(X){ return( -sin(2*X[,1]) ) } 
g2 = function(X){ return( cos(3*X[,2])  ) } 
int = function(X){
  x1 = X[,1]
  x2 = X[,2]
  return( 2*cos(x1)*sin(2*x2)  ) 
}


m_full_why = function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return( -sin(2*x1) + cos(3*x2) + 2*cos(x1)*sin(2*x2) )
}


true_model_list = list()
true_model_list[[1]] = m_x1
true_model_list[[2]] = m_x2
true_model_list[[3]] = m_full_why 


l = -2; u = 2
l_int = l; u_int = u
d <<- 2

samp = vector(mode="list")

#generate B bootstrap samples
ci_points = seq(-2, 2, length.out=100)
bs_shap=matrix(0, ncol=2, nrow=B)
bs_shap_int=matrix(0, ncol=2, nrow=B)

for ( i in 1:B){
  
model_list_bs = vector(mode="list")
samp = dat[sample(nrow(dat), N, replace=TRUE), ]
model.np1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = samp, bws=b1)
model.np2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = samp, bws=b2) 
model.np12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = samp, bws=b12 ) 
model_list_bs[[1]] = model.np1 
model_list_bs[[2]] = model.np2 
model_list_bs[[3]] = model.np12 
bs_shap[i,1] = shapley_vec(j=1, c(ci_points[j],-0.5), model_list = model_list_bs) #component-based x1
bs_shap[i,2] = shapley_vec(j=2, c(ci_points[j],-0.5), model_list = model_list_bs) #component-based x2

bs_shap_int[i,1] = shapley_int(j=1, c(ci_points[j],-0.5), model_list = model_list_bs, 
                               model_list_int = model_list_int) #integration-based x1
bs_shap_int[i,2] = shapley_int(j=2, c(ci_points[j],-0.5), model_list = model_list_bs, 
                               model_list_int = model_list_int) #integration-based x2

}

alpha = 0.05
qs1 = quantile(bs_shap[,1], probs=c(alpha/2, 1 - alpha/2) )
qs2 = quantile(bs_shap[,2], probs=c(alpha/2, 1 - alpha/2) )

qs1_int = quantile(bs_shap_int[,1], probs=c(alpha/2, 1 - alpha/2) )
qs2_int = quantile(bs_shap_int[,2], probs=c(alpha/2, 1 - alpha/2) )


return(c(qs1, qs2, qs1_int, qs2_int))
}


library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)


source("functions.R")
N=2000
cov=0
sigma_sim<<-matrix(c(4, cov,
                     cov, 4), nrow=2, ncol=2)


source("integral_population_2d_int.R")
source("integral_estimation_2d.R") 
source("shapley_int.R")

g1 = function(X){ return( -sin(2*X[,1]) ) } 
g2 = function(X){ return( cos(3*X[,2])  ) } 
int = function(X){
  x1 = X[,1]
  x2 = X[,2]
  return( 2*cos(x1)*sin(2*x2)  ) 
}




m_full_why = function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return( -sin(2*x1) + cos(3*x2) + 2*cos(x1)*sin(2*x2) )
}


true_model_list = list()
true_model_list[[1]] = m_x1
true_model_list[[2]] = m_x2
true_model_list[[3]] = m_full_why 


l = -2; u = 2 
l_int = l; u_int = u
d <<- 2

Y=dat[,1]
#initial estimation
model.np.init1 <<- npreg(reformulate("X1", "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
model.np.init2 <<- npreg(reformulate("X2", "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
model.np.init12 <<- npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
#shapley estimator
model_list_init <<- vector(mode="list")
model_list_init[[1]] = model.np.init1 
model_list_init[[2]] = model.np.init2
model_list_init[[3]] = model.np.init12 
#take bandwidths
b1 <<- model.np.init1$bw
b2 <<- model.np.init2$bw
b12 <<- model.np.init12$bw


library(parallel)
qs = mclapply(1:100, bootstrap, mc.cores=40)

points=matrix(unlist(qs), byrow = FALSE, ncol=100)



#Plots:
x1_grid = seq(-2, 2, length.out=1000)  
x2 = rep(-0.5, 1000)
coordinates = rbind(x1_grid, x2) 

shap_eval1 = matrix(0, nrow=ncol(coordinates), ncol=1)
shap_eval2 = matrix(0, nrow=ncol(coordinates), ncol=1)

for (i in 1:ncol(coordinates)){
  print(i)
  grid_col=as.numeric(coordinates[,i])
  shap_eval1[i] = shapley_popul(j=1, grid_col, model_list = model_list_init) 
  shap_eval2[i] = shapley_popul(j=2, grid_col, model_list = model_list_init) 
  
}

#comp estimator
shap_eval_est1=shapley_vec(j=1, coordinates, model_list = model_list_init) 
shap_eval_est2=shapley_vec(j=2, coordinates, model_list = model_list_init)

#integration estimator
shap_eval_est1 = matrix(0, nrow=ncol(coordinates), ncol=1)
shap_eval_est2 = matrix(0, nrow=ncol(coordinates), ncol=1)

for (i in 1:ncol(coordinates)){
  print(i)
  grid_col=as.numeric(coordinates[,i])
  shap_eval_est1[i] = shapley_int(j=1, grid_col, model_list = model_list_init, model_list_int = model_list_int) 
  shap_eval_est2[i] = shapley_int(j=2, grid_col, model_list = model_list_init, model_list_int = model_list_int) 
  
}

x1_grid_ex = seq(-2, 2, length.out=100)
plot(y=shap_eval1, x=x1_grid, type="l", xlab="x1", ylab="")
points(y=shap_eval_est1, x=x1_grid, col="blue", type="l")
points(y = points[5,], x = x1_grid_ex, col="red", cex=0.3)
points(y = points[6,], x = x1_grid_ex, col="red", cex=0.3)


plot(y=shap_eval2, x=x1_grid, type="l", xlab="x1", ylab="", ylim=c(-1,0.5))
points(y=shap_eval_est2, x=x1_grid, col="blue", type="l")
points(y = points[7,], x = x1_grid_ex, col="red", cex=0.3)
points(y = points[8,], x = x1_grid_ex, col="red", cex=0.3)







