m_full= function(x1, x2, x3){
  return(-sin(2*x1) + cos(3*x2) + 2*cos(x1)*sin(2*x2) + 0.5*x3)
} 

m_full_why <<- function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return(-sin(2*x1) + cos(3*x2) + 2*cos(x1)*sin(2*x2) + 0.5*x3)
}


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


norm2_vec = function(dep, cond){
  x = dep # dependent variable, 2D vector
  Y_2 = cond # what you condition on, skalar
  c_mu_2 = inv_yy_2 %*% Y_2
  end = ncol(x)  
  
  return(sapply(1:end, 
         function(a){
          pre_mult2*exp(-0.5*t(x[,a] - c_mu_2) %*% inv %*% (x[,a] - c_mu_2))
         })) # returns normal draws in vector

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



# indep and dep 
m_x3_temp=function(x_out, x3){

  x1=x_out[1,]
  x2=x_out[2,]

    res=m_full(x1, x2, x3)*norm2_vec(dep=x_out, cond=x3) # vectorized normal
  return(res)
}

m_x3=function(X){
  x3=as.numeric(X[3])
  return(cubintegrate(f = m_x3_temp, lower = rep(-5, 2), upper = rep(5, 2), method = "cuhre",
               relTol = 3e-1, x3=x3, nVec = 128L)$integral)
}# m(x3)






m_x2_temp=function(x_out, x2){
  x1=x_out[1,]
  x3=x_out[2,]
  
  res=m_full(x1, x2, x3)*norm2_vec(dep=x_out, cond=x2)
  return(res)
}

m_x1_temp=function(x_out, x1){
  x2=x_out[1,]
  x3=x_out[2,]
  
  res=m_full(x1, x2, x3)*norm2_vec(dep=x_out, cond=x1)
  return(res)
}


m_x2=function(X){
  x2=as.numeric(X[2])
  
  return(cubintegrate(f = m_x2_temp, lower = rep(-5, 2), upper = rep(5, 2), method = "cuhre",
               relTol = 3e-1, x2=x2, nVec = 128L)$integral)
}# m(x2)
m_x1=function(X){
  x1=as.numeric(X[1])
  return(cubintegrate(f = m_x1_temp, lower = rep(-5, 2), upper = rep(5, 2), method = "cuhre",
                      relTol = 3e-1, x1=x1, nVec = 128L)$integral)
}# m(x1)






m_x1_x3_temp=function(x_out, x1, x3){
  x2=x_out
  
  res = m_full(x1, x2, x3)*norm1_vec(dep=x2, cond=c(x1,x3))
  
  return(res)
} # m(x1,x3)


m_x1_x3=function(X){
  x1=as.numeric(X[1])
  x3=as.numeric(X[3])
  return(cubintegrate(f = m_x1_x3_temp, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x1=x1, x3=x3, nVec = 128L)$integral)
  
}# m(x1,x3)







m_x2_x3_temp=function(x_out, x2, x3){
  x1=x_out

  res = m_full(x1, x2, x3)*norm1_vec(dep=x1, cond=c(x2,x3))
  
  return(res)
} # m(x2,x3)



m_x1_x2_temp=function(x_out, x1, x2){
  x3=x_out
  
  res = m_full(x1, x2, x3)*norm1_vec(dep=x3, cond=c(x1,x2))
  
  return(res)
} # m(x1,x2)



m_x1_x2=function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  return(cubintegrate(f = m_x1_x2_temp, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x1=x1, x2=x2, nVec = 128L)$integral) 
  }# m(x1,x2)





m_x2_x3=function(X){
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return(cubintegrate(f = m_x2_x3_temp, lower = rep(-5, 1), upper = rep(5, 1), method = "cuhre",
                      relTol = 3e-1, x2=x2, x3=x3, nVec = 128L)$integral)
  }# m(x2,x3)




true_model_list = list()
true_model_list[[1]] = m_x1
true_model_list[[2]] = m_x2
true_model_list[[3]] = m_x3 
true_model_list[[4]] = m_x1_x2
true_model_list[[5]] = m_x1_x3
true_model_list[[6]] = m_x2_x3 
true_model_list[[7]] = m_full_why 


