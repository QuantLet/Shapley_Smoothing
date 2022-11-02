

#cov=0.3
#sigma=matrix(c(4, cov, cov,
#               cov, 4, cov,
#               cov, cov, 4), nrow=3, ncol=3)


# Can be used for x1/x2,x3 and x2/x1,x3 and x3/x1,x2
sigma_xy_1 = t(sigma_sim[1, 2:3])
sigma_yy_1 = sigma_sim[2:3, 2:3]
sigma_yx_1 = t(sigma_xy_1)
sigma_xx_1 = sigma_sim[1,1]
#conditional variance
inv_yy_1 = sigma_xy_1 %*% solve(sigma_yy_1)
c_var_1 = sigma_xx_1 - sigma_xy_1 %*% solve(sigma_yy_1) %*% sigma_yx_1
pre_mult = 1/(sqrt(c_var_1*2*pi))
sq = sqrt(c_var_1)

norm1 = function(dep, cond){
  
  x = dep # dependent variable, skalar
  Y_1 = cond # what you condition on, 2D vector
  c_mu_1 = inv_yy_1 %*% Y_1
  #cond normal density:
  return(pre_mult* exp(  -0.5*((x - c_mu_1)/ sq  )^2 ))
}


#Can be used for x1,x2/x3 and x2,x3/x1 and x1,x3/x2
sigma_xy_2 = matrix(sigma_sim[1:2, 3])
sigma_yy_2 = sigma_sim[3, 3]
sigma_yx_2 = t(sigma_xy_2)
sigma_xx_2 = sigma_sim[1:2, 1:2]
#conditional variance
inv_yy_2 = sigma_xy_2 %*% solve(sigma_yy_2)
c_var_2 = sigma_xx_2 - sigma_xy_2 %*% as.numeric(solve(sigma_yy_2)) %*% sigma_yx_2
pre_mult2 = (det(2*pi*sigma_xx_2))^(-0.5)
inv = solve(sigma_xx_2)

norm2 = function(dep, cond){
  x = dep # dependent variable, 2D vector
  Y_2 = cond # what you condition on, skalar
  c_mu_2 = inv_yy_2 %*% Y_2
  
  #cond multivariate normal density:
  return(pre_mult2*exp(-0.5*t( x - c_mu_2  ) %*% inv %*% (x - c_mu_2) ))
}





norm1_vec = function(dep, cond){
  
  x = dep # dependent variable, skalar
  Y_1 = cond # what you condition on, 2D vector
  c_mu_1 = inv_yy_1 %*% Y_1
  #cond normal density:
  
  return(sapply(1:length(x), function(a){
    
    pre_mult*exp(-0.5*((x[a] - c_mu_1)/sq)^2)
    
  } ))
  #return(pre_mult* exp(  -0.5*((x - c_mu_1)/ sq  )^2 ))
}

#a=seq(-4,4, length.out=8)
#for (i in 1:length(a)){
#  print(norm1(dep=a[i], cond=c(2,2)))
#  print(dcmvnorm(x=a[i], mean=rep(0, 3), sigma=sigma, log=FALSE, dep=c(3), given=c(1,2), X=c(2,2))) # x1,x2/x3
#  print(dnorm(a[i], 0, 2))

#}

#b=cbind(seq(-4,4, length.out=8), seq(-4,4, length.out=8))
#for (i in 1:length(a)){
#  print(norm2(dep=b[i,], cond=0 ))
#  print(dcmvnorm(x=b[i,], mean=rep(0, 3), sigma=sigma, log=FALSE, dep=c(1,2), given=c(3), X=c(0)) # x1,x2/x3
#  )
#}

norm2_vec = function(dep, cond){
  x = dep # dependent variable, 2D vector
  Y_2 = cond # what you condition on, skalar
  c_mu_2 = inv_yy_2 %*% Y_2
  #end=ncol(as.matrix(x))
  end = ncol(x)
  
  return(sapply(1:end, 
                function(a){
                  pre_mult2*exp(-0.5*t(x[,a] - c_mu_2) %*% inv %*% (x[,a] - c_mu_2))
                })) # returns normal draws in vector
  #cond multivariate normal density:
  #  return(pre_mult2*exp(-0.5*t( x - c_mu_2  ) %*% inv %*% (x - c_mu_2) ))
}


# indep and dep 
m_x3_temp_int=function(x_out, x3){
  x1=x_out[1,]
  x2=x_out[2,]
  
  #res=m_full(x1, x2, x3)*dnorm(x1,0,2)*dnorm(x2,0,2)
  #res= m_full(x1, x2, x3)*dcmvnorm(x=c(x1, x2), mean=rep(0, d), 
  #                                sigma=cmat, log=FALSE, dep=c(1,2), given=c(3), X=c(x3)) # x1,x2/x3
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
  m_full_hat = predict(model_list[[7]], newdata = x_eval_vec) 
  #dep = rbind(x1,x2)
  res = m_full_hat*norm2_vec(dep=x_out, cond=x3)
  return(res)
}


m_x3_est=function(X){
  x3=as.numeric(X[3])
  #  return(hcubature(f = m_x3_temp, rep(-5, 2), rep(5, 2), tol=1e-1, x3=x3)$integral)
  return(cubintegrate(f = m_x3_temp_int, lower = rep(-5, 2), upper = rep(5, 2), method = "cuhre",
                      relTol = 3e-1, x3=x3, nVec = 128L)$integral)
}# m(x3)



 # X = data.frame(X1=seq(-2,2,length.out=500), X2=seq(-2,2,length.out=500), X3 = seq(-2,2,length.out=500) )
#  y2 = g2(X)
#  y = sapply( 1:500, function(a){ 
#    print(a)
#    return(m_x2_est(X[a,]))} )


#  plot(y2, x=seq(-2,2,length.out=500), type="l")
#  lines(y, x=seq(-2,2,length.out=500), type="l", col="blue")



m_x2_temp_int=function(x_out, x2){
  x1=x_out[1,]
  x3=x_out[2,]
  
  #res=m_full(x1, x2, x3)*dnorm(x1,0,2)*dnorm(x3,0,2)
  
  #res = m_full(x1, x2, x3)*dcmvnorm(x=c(x1, x3), mean=rep(0, d), 
  #                                sigma=cmat, log=FALSE, dep=c(1,3), given=c(2), X=c(0)) # x1,x3/x2
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
  m_full_hat = predict(model_list[[7]], newdata = x_eval_vec) 
  #dep=rbind(x1,x3)
  res = m_full_hat*norm2_vec(dep=x_out, cond=x2)
  return(res)
}

m_x1_temp_int=function(x_out, x1){
  x2=x_out[1,]
  x3=x_out[2,]
  
  #res=m_full(x1, x2, x3)*dnorm(x2,0,2)*dnorm(x3,0,2)
  
  #res = m_full(x1, x2, x3)*dcmvnorm(x=c(x2, x3), mean=rep(0, d), 
  #                                 sigma=cmat, log=FALSE, dep=c(2,3), given=c(1), X=c(0)) # x2,x3/x1
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
  m_full_hat = predict(model_list[[7]], newdata = x_eval_vec) 
  
  #dep = rbind(x2,x3)
  res = m_full_hat*norm2_vec(dep=x_out, cond=x1)
  return(res)
}



m_x1_x3_temp_int=function(x_out, x1, x3){
  x2=x_out
  
  #res=m_full(x1, x2, x3)*dnorm(x2,0,2)
  #res = m_full(x1, x2, x3)*dcmvnorm(x=c(x2), mean=rep(0, d), 
  #                                  sigma=cmat, log=FALSE, dep=c(2), given=c(1,3), X=c(0,0)) # x2/x1,x3
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
  m_full_hat = predict(model_list[[7]], newdata = x_eval_vec) 
  res = m_full_hat*norm1_vec(dep=x2, cond=c(x1,x3))
  
  return(res)
} # m(x1,x3)

m_x2_x3_temp_int=function(x_out, x2, x3){
  x1=x_out
  
  #res=m_full(x1, x2, x3)*dnorm(x1,0,2)
  # res = m_full(x1, x2, x3)*dcmvnorm(x=c(x1), mean=rep(0, d), 
  #                                    sigma=cmat, log=FALSE, dep=c(1), given=c(2,3), X=c(0,0)) # x1/x2,x3
  
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
  m_full_hat = predict(model_list[[7]], newdata = x_eval_vec) 
  res = m_full_hat*norm1_vec(dep=x1, cond=c(x2,x3))
  
  return(res)
} # m(x2,x3)



m_x1_x2_temp_int=function(x_out, x1, x2){
  x3=x_out 
  
  #res=m_full(x1, x2, x3)*dnorm(x3,0,2)
  
  #res = m_full(x1, x2, x3)*dcmvnorm(x=c(x3), mean=rep(0, d), 
  #                                  sigma=cmat, log=FALSE, dep=c(3), given=c(1,2), X=c(0,0)) # x3/x1,x2
  
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
  m_full_hat = predict(model_list[[7]], newdata = x_eval_vec) 
  
  res = m_full_hat*norm1_vec(dep=x3, cond=c(x1,x2))
  
  return(res)
} # m(x1,x2)



m_x1_x2_est=function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
#  return(hcubature(f= m_x1_x2_temp, rep(-5, 1), rep(5, 1), tol=1e-1, x1=x1, x2=x2)$integral)
  return(cubintegrate(f = m_x1_x2_temp_int, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x1=x1, x2=x2, nVec = 128L)$integral) 
}# m(x1,x2)


m_x1_x3_est=function(X){
  x1=as.numeric(X[1])
  x3=as.numeric(X[3])
 # return(hcubature(f= m_x1_x3_temp, rep(-5, 1), rep(5, 1), tol=1e-1, x1=x1, x3=x3)$integral)
  return(cubintegrate(f = m_x1_x3_temp_int, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x1=x1, x3=x3, nVec = 128L)$integral)
}# m(x1,x3)

m_x2_x3_est=function(X){
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
 # return(hcubature(f= m_x2_x3_temp, rep(-5, 1), rep(5, 1), tol=1e-1, x2=x2, x3=x3)$integral)
  return(cubintegrate(f = m_x2_x3_temp_int, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x2=x2, x3=x3, nVec = 128L)$integral)
}# m(x2,x3)





m_x2_est=function(X){
  x2=as.numeric(X[2])
  return(cubintegrate(f = m_x2_temp_int, lower = rep(-5, 2), upper = rep(5, 2), method = "cuhre",
                      relTol = 3e-1, x2=x2, nVec = 128L)$integral)
  # return(hcubature(f= m_x2_temp, rep(-5, 2), rep(5, 2), tol=1e-1, x2=x2)$integral)
}# m(x2)
m_x1_est=function(X){
  x1=as.numeric(X[1])
  return(cubintegrate(f = m_x1_temp_int, lower = rep(-5, 2), upper = rep(5, 2), method = "cuhre",
                      relTol = 3e-1, x1=x1, nVec = 128L)$integral)
#  return(hcubature(f= m_x1_temp, rep(-5, 2), rep(5, 2), tol=1e-1, x1=x1)$integral)
}# m(x1)


m_full_hat = function(X){
 #x1 = as.numeric(X[1])
 #x2 = as.numeric(X[2])
 #x3 = as.numeric(X[3])
 #x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2), X3 = c(x3) )  
   return(predict(model_list[[7]], newdata = X) ) # X is a df
}
  


model_list_int = list(m_x1_est, m_x2_est, m_x3_est, m_x1_x2_est, m_x1_x3_est, m_x2_x3_est, m_full_hat)


