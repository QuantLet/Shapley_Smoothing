shapley_vec = function(j, x_eval, model_list){
#print(x_eval)
  x_eval = t(x_eval)
  shap = rep(0, nrow(x_eval)) 
  x_eval_vec <- data.frame(X1 = c(x_eval[,1]), X2 = c(x_eval[,2]) )  
  
  #new begin
  c = mean(Y) # empty set
  shap =  -(1/d)*((nchoosek(d-1, 0))^(-1))*rep(c, nrow(x_eval) )
  #new end
  
  for(k in 1:length(model_list) ){
    pred = predict(model_list[[k]], newdata = x_eval_vec[, c(model_list[[k]]$xnames), drop = FALSE ]) 
    shap = shap + weight(j,k,model_list)*pred
  }
  
  return(as.matrix(shap))
}