library(parallel)
  
bs=function(m){
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
    return(-sin(2*x1) + cos(3*x2) + 2*cos(x1)*sin(2*x2))
  }
  
  
  true_model_list = list()
  true_model_list[[1]] = m_x1
  true_model_list[[2]] = m_x2
  true_model_list[[3]] = m_full_why 
  
  l = -2; u = 2; N=4000
  l_int = l; u_int = u
  
  d <<- 2
  
X<<-data.frame(mvrnorm(n=N, mu=c(0,0), Sigma=sigma_sim[1:2,1:2]))
#DGP

Y <<- g1(X) + g2(X) + int(X) + rnorm(n=nrow(X))
dat = data.frame(Y,X)
samp = vector(mode="list")

#generate B bootstrap samples
B=1000
bs_shap=numeric()
bs_shap2=numeric()

bs_shap_int=numeric()
bs_shap_int2=numeric()

#initial estimation
model.np.init1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
model.np.init2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")
model.np.init12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = dat, nmulti=1, bwmethod="cv.aic")

#shapley estimator
model_list_init = vector(mode="list")
model_list_init[[1]] = model.np.init1 
model_list_init[[2]] = model.np.init2
model_list_init[[3]] = model.np.init12 

#shap_estimated = shapley_vec(j=1, c(0,0), model_list = model_list_init) 
#shap_estimated2 = shapley_vec(j=2, c(0,0), model_list = model_list_init) 

#take bandwidths
b1 = model.np.init1$bw
b2 = model.np.init2$bw
b12 = model.np.init12$bw

for ( i in 1:B){
model_list_bs = vector(mode="list")
samp = dat[sample(nrow(dat), N, replace=TRUE), ]


model.np1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = samp, bws=b1)
model.np2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = samp, bws=b2) 
model.np12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = samp, bws=b12 ) 

model_list_bs[[1]] = model.np1 
model_list_bs[[2]] = model.np2 
model_list_bs[[3]] = model.np12 

#bs_shap[i] = shapley_vec(j=1, c(0,0), model_list = model_list_bs) #component-based x1
#bs_shap2[i] = shapley_vec(j=2, c(0,0), model_list = model_list_bs) #component-based x2

bs_shap_int[i] = shapley_int(j=1, c(-0.5,-0.5), model_list = model_list_bs, 
                               model_list_int = model_list_int) #integration-based x1
bs_shap_int2[i] = shapley_int(j=2, c(-0.5,-0.5), model_list = model_list_bs, 
                               model_list_int = model_list_int) #integration-based x2

}

#quantiles
# take alpha/2 and 1 - alpha/2 quantile
alpha1 = 0.15
qs1 = quantile(bs_shap_int, prob=c(alpha1/2, 1 - alpha1/2))
qs1_2 = quantile(bs_shap_int2, prob=c(alpha1/2, 1 - alpha1/2))

alpha2 = 0.10
qs2 = quantile(bs_shap_int, prob=c(alpha2/2, 1 - alpha2/2))
qs2_2 = quantile(bs_shap_int2, prob=c(alpha2/2, 1 - alpha2/2))

alpha3 = 0.05
qs3 = quantile(bs_shap_int, prob=c(alpha3/2, 1 - alpha3/2))
qs3_2 = quantile(bs_shap_int2, prob=c(alpha3/2, 1 - alpha3/2))

return(c(qs1, qs2, qs3, qs1_2, qs2_2, qs3_2))
}



M=1000
results_list = mclapply(1:M, bs, mc.cores=48)

results = matrix(unlist(results_list), byrow = FALSE, ncol=1000)






