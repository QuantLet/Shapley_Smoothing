shapley_popul = function(j, x_eval){
  shap = 0  
  
  shap_res=sapply(1:7, function(k){
    shap = weight(j,k)*true_model_list[[k]](x_eval)
    
  }
  )
  
  return(sum(shap_res))
}




