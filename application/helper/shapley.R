shapley = function(j, x_eval){
  shap = 0
  x_eval_vec <- data.frame(X1 = c(x_eval[1]), X2 = c(x_eval[2]), X3 = c(x_eval[3]) )  
  c = mean(Y) 
  shap =  -(1/d)*((nchoosek(d-1, 0))^(-1))*c
  for(k in 1:length(model_list) ){
    pred = predict(model_list[[k]], newdata = x_eval_vec[, c(model_list[[k]]$xnames), drop = FALSE ]) 
    shap = shap + weight(j,k)*pred
  }
  
  return(as.matrix(shap))
}


