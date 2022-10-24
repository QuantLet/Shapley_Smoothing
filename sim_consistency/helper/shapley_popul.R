shapley_popul = function(j, x_eval){
  #shap = matrix(0, nrow=nrow(x_eval)) #if matrix evaluation points
  shap = 0  
  #x_eval <- data.frame(X1 = c(x_eval[1]), X2 = c(x_eval[2]), X3 = c(x_eval[3]))  
  
  #length(true_model_list) is 7
 # for(k in 1:7 ){
  #  pred = true_model_list[[k]](as.numeric(x_eval))
  #  shap = shap + weight(j,k)*pred
  #}
  
  shap_res=sapply(1:7, function(k){
    shap = weight(j,k)*true_model_list[[k]](x_eval)
    
  }
  )
  
  return(sum(shap_res))
}




