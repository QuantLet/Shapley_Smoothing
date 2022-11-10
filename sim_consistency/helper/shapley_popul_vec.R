shapley_popul_vec = function(j, x_eval){
  x_eval = t(x_eval)
  shap = rep(0, nrow(x_eval)) 
  x_eval_vec <- data.frame(X1 = c(x_eval[,1]), X2 = c(x_eval[,2]), X3 = c(x_eval[,3]))  
  
  for(k in 1:length(true_model_list) ){
    
    pred = true_model_list[[k]](x_eval_vec)
    
    shap = shap + weight(j,k)*pred
  }
  
  return(as.matrix(shap))
}