install.packages("remotes")
remotes::install_github("reyesem/IntroAnalysis")
  
bs=function(m){
  library(np)
  library(pracma)
  library(cubature)
  library(simstudy)
  library(MASS)
  library(IntroAnalysis)

  
  source("functions.R")
  
  cov=0
  sigma_sim<<-matrix(c(4, cov,
                     cov, 4), nrow=2, ncol=2)
  
  
  source("integral_population_2d.R")
  source("integral_estimation_2d.R") 
  source("shapley_int.R")

  
  g1 = function(X){ return( -sin(2*X[,1]) ) } #E(g1)=0. 
  g2 = function(X){ return( 0.1*X[,2]  ) } #E(g2)=0 
  int = function(X){
    x1 = X[,1]
    x2 = X[,2]
    return( 2*cos(x1)*sin(x2)  ) 
  }
  
  
  m_full_why = function(X){
    x1=as.numeric(X[1])
    x2=as.numeric(X[2])
    return(-sin(2*x1) + 0.1*x2 + 2*cos(x1)*sin(x2))
  }
  
  
  true_model_list = list()
  true_model_list[[1]] = m_x1
  true_model_list[[2]] = m_x2
  true_model_list[[3]] = m_full_why 
  
  
  l = -2; u = 2; N=500; 
  l_int = l; u_int = u
  
  d <<- 2
  

X<<-data.frame(mvrnorm(n=N, mu=c(0,0), Sigma=sigma_sim[1:2,1:2]))

#DGP
Y <<- g1(X) + g2(X) + int(X) + rnorm(nrow(X))
dat = data.frame(Y,X)
samp = vector(mode="list")

#generate B bootstrap samples
B=1000
pred_bs=numeric()
bs_shap = numeric()
bs = numeric()
bs2 = numeric()


#initial estimation
model.np.init1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
model.np.init2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
model.np.init12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")

#shapley estimator
model_list_init = vector(mode="list")
model_list_init[[1]] = model.np.init1 
model_list_init[[2]] = model.np.init2
model_list_init[[3]] = model.np.init12 

#component-based x1
shap_estimated = shapley_vec(j=1, c(0, 0), model_list = model_list_init) 

#component-based x2
shap_estimated2 = shapley_vec(j=2, c(0, 0), model_list = model_list_init) 

#take bandwidths
b1 = model.np.init1$bw
b2 = model.np.init2$bw
b12 = model.np.init12$bw

#modify bandwidth of full model:
g12 = b12 * log(log(N))*4
g1 = b1 * log(log(N))*4
g2 = b2 * log(log(N))*4


mod_g = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = dat, bws=g12)
pred_g = predict(mod_g, newdata = X)

mod_g1 = npreg(reformulate(c("X1"), "Y"), regtype = "ll", data = dat, bws=g1) 
pred1 = predict(mod_g1, newdata = X)

mod_g2 = npreg(reformulate(c("X2"), "Y"), regtype = "ll", data = dat, bws=g2) 
pred2 = predict(mod_g2, newdata = X)

df_pred = data.frame(X1 = c(0), X2 = c(0))
mod_h = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = dat, bws=b12) 


#Step 1: Construct residuals
eps = Y - predict(model_list_init[[3]], newdata = X) 
eps1 = Y - predict(model_list_init[[1]], newdata = X) 
eps2 = Y - predict(model_list_init[[2]], newdata = X) 

pred_g_point = predict(mod_g, newdata =  df_pred)
pred_g1_point = predict(mod_g1, newdata =  df_pred)
pred_g2_point = predict(mod_g2, newdata =  df_pred)


for (i in 1:B){
  
  
#step 2: Multiply with rv V
random = rmammen(n=N, construct= c("two-point mass"))

#step 3: construct y_b
Y_b = pred_g + eps*random
Y_b1 = pred1 + eps1*random
Y_b2 = pred2 + eps2*random

model.np1 = npreg(reformulate(c("X1"), "Y_b1"), regtype = "ll", data = data.frame(Y_b1, X), bws=b1) 
model.np2 = npreg(reformulate(c("X2"), "Y_b2"), regtype = "ll", data = data.frame(Y_b2, X), bws=b2) 
model.np12 = npreg(reformulate(c("X1", "X2"), "Y_b"), regtype = "ll", data = data.frame(Y_b, X), bws=b12) 

pred_bs12 = predict(model.np12, newdata = df_pred)
pred_bs1 = predict(model.np1, newdata = df_pred)
pred_bs2 = predict(model.np2, newdata = df_pred)


w11 = weight(j=1, k=1, model_list = model_list_init)
w12 = weight(j=1, k=2, model_list = model_list_init)
w_empty = -(1/d)*((nchoosek(d-1, 0))^(-1))

bs[i] = (1/d)*(pred_bs12 - pred_g_point) + w11*(pred_bs1 - pred_g1_point) + w12*(pred_bs2 - pred_g2_point) +
  w_empty*(mean(Y_b) - mean(Y))

bs2[i] = (1/d)*(pred_bs12 - pred_g_point) + w12*(pred_bs1 - pred_g1_point) + w11*(pred_bs2 - pred_g2_point) +
  w_empty*(mean(Y_b) - mean(Y))
  

}


qs1 = quantile(bs, prob=c(0.05/2, 1 - 0.05/2))
qs2 = quantile(bs, prob=c(0.1/2, 1 - 0.1/2))
qs3 = quantile(bs, prob=c(0.15/2, 1 - 0.15/2))

qs1x2 = quantile(bs2, prob=c(0.05/2, 1 - 0.05/2))
qs2x2 = quantile(bs2, prob=c(0.1/2, 1 - 0.1/2))
qs3x2 = quantile(bs2, prob=c(0.15/2, 1 - 0.15/2))


 CI_lower1 = shap_estimated + qs1[1]
 CI_upper1 = shap_estimated + qs1[2]
 qs1_res = c(CI_lower1, CI_upper1)
 
 CI_lower2 = shap_estimated + qs2[1]
 CI_upper2 = shap_estimated + qs2[2]
 qs2_res = c(CI_lower2, CI_upper2)
 
 CI_lower3 = shap_estimated + qs3[1]
 CI_upper3 = shap_estimated + qs3[2]
 qs3_res = c(CI_lower3, CI_upper3)

 
 CI2_lower1 = shap_estimated2 + qs1x2[1]
 CI2_upper1 = shap_estimated2 + qs1x2[2]
 qs1_res2 = c(CI2_lower1, CI2_upper1)
 
 CI2_lower2 = shap_estimated2 + qs2x2[1]
 CI2_upper2 = shap_estimated2 + qs2x2[2]
 qs2_res2 = c(CI2_lower2, CI2_upper2)
 
 CI2_lower3 = shap_estimated2 + qs3x2[1]
 CI2_upper3 = shap_estimated2 + qs3x2[2]
 qs3_res2 = c(CI2_lower3, CI2_upper3)


result = c(qs3_res, qs2_res, qs1_res, qs3_res2, qs2_res2, qs1_res2)

return(result)
}

library(parallel)
M=1000
results_list = mclapply(1:M, bs, mc.cores=48)
results = matrix(unlist(results_list), byrow = FALSE, ncol=M)
qs1 = t(results)

shap_true = -0.002487695 # for x1

sum(shap_true >= qs1[,1] & shap_true <= qs1[,2])/M 
sum(shap_true >= qs1[,3] & shap_true <= qs1[,4])/M 
sum(shap_true >= qs1[,5] & shap_true <= qs1[,6])/M 

shap_true = -0.002487695# for x2

sum(shap_true >= qs1[,7] & shap_true <= qs1[,8])/M 
sum(shap_true >= qs1[,9] & shap_true <= qs1[,10])/M 
sum(shap_true >= qs1[,11] & shap_true <= qs1[,12])/M 


shap_true = -0.7727975 # for x2
sum(shap_true >= qs1[,13] & shap_true <= qs1[,14])/M 
sum(shap_true >= qs1[,15] & shap_true <= qs1[,16])/M 
sum(shap_true >= qs1[,17] & shap_true <= qs1[,18])/M 
sum(shap_true >= qs1[,19] & shap_true <= qs1[,20])/M 
sum(shap_true >= qs1[,21] & shap_true <= qs1[,22])/M 
sum(shap_true >= qs1[,23] & shap_true <= qs1[,24])/M 


