shapley_int = function(j, x_eval, model_list, model_list_int){
  shap = 0
  model_list <<- model_list # global 
  
  x_eval_vec <- data.frame(X1 = c(x_eval[1]), X2 = c(x_eval[2]))  
  
  c = mean(Y) # empty set
  shap =  -(1/d)*((nchoosek(d-1, 0))^(-1))*c
  for(k in 1:length(model_list) ){
    pred = model_list_int[[k]](x_eval_vec)
    shap = shap + weight(j,k,model_list)*pred
  }
  

    return(shap)
}


#hcubature(f=shapley_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=2)
















