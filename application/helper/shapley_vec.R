shapley_vec = function(j, x_eval, alt, model_list){

  x_eval = t(x_eval)
  shap = rep(0, nrow(x_eval)) 
  x_eval_vec <- data.frame(X1 = c(x_eval[,1]), X2 = c(x_eval[,2]), X3 = c(x_eval[,3]) )  
  c = mean(Y) 
  shap =  -(1/d)*((nchoosek(d-1, 0))^(-1))*rep(c, nrow(x_eval) )

  if(alt==FALSE){
  for(k in 1:length(model_list) ){
    names_vars = model_list[[k]]$xnames 
    pred = predict(model_list[[k]], newdata = x_eval_vec[, c(names_vars), drop = FALSE ]) 
    shap = shap + weight(j, k, names_vars)*pred
  }
  } else{
    for(k in 1:length(model_list) ){
      names_vars = row.names(model_list[[k]]$importance)
      pred = predict(model_list[[k]], newdata = x_eval_vec[, c(names_vars), drop = FALSE ]) 
      shap = shap + weight(j,k, names_vars)*pred
    }
}
  
  
  return(as.matrix(shap))
}