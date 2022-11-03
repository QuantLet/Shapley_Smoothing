


# Can be used for x1/x2,x3 and x2/x1,x3 and x3/x1,x2
sigma_xy_1 = t(sigma_sim[1, 2])
sigma_yy_1 = sigma_sim[2, 2]
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

norm1_vec = function(dep, cond){
  
  x = dep # dependent variable, skalar
  Y_1 = cond # what you condition on, 2D vector
  c_mu_1 = inv_yy_1 %*% Y_1
  #cond normal density:
  
  return(sapply(1:length(x), function(a){
    
    pre_mult*exp(-0.5*((x[a] - c_mu_1)/sq)^2)
    
  } ))
}




#x1/x2
m_x2_temp_int=function(x_out, x2){
  x1=x_out
  
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2) )  
  m_full_hat = predict(model_list[[3]], newdata = x_eval_vec) 
  
  
  res = m_full_hat*norm1_vec(dep=x1, cond=c(x2))
  
  return(res)
} # m(x2)


m_x2_est=function(X){
  x2=as.numeric(X[2])
  return(cubintegrate(f = m_x2_temp_int, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x2=x2, nVec = 128L)$integral)
  
}# m(x2)

#x2/x1
m_x1_temp_int=function(x_out, x1){
  x2=x_out
  
  x_eval_vec = data.frame(X1 = c(x1), X2 = c(x2) )  
  m_full_hat = predict(model_list[[3]], newdata = x_eval_vec) 
  res = m_full_hat*norm1_vec(dep=x2, cond=c(x1))
  
  return(res)
} # m(x1)


m_x1_est=function(X){
  x1=as.numeric(X[1])
  return(cubintegrate(f = m_x1_temp_int, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x1=x1, nVec = 128L)$integral)
  
}# m(x1)









m_full_hat = function(X){

   return(predict(model_list[[3]], newdata = X) ) # X is a df
}
  


model_list_int = list(m_x1_est, m_x2_est, m_full_hat)



