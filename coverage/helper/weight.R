weight = function(j,k, model_list){
  
  indicator = as.integer(sum(as.integer(paste("X",j, sep="") == c(model_list[[k]]$xnames)))>0) # j=variable, k= subset s
  sign = ifelse(indicator > 0, 1, -1) # sign fct
  
  card_s = model_list[[k]]$ndim
  
  return( sign*(1/d)*( nchoosek(d-1, card_s-indicator) )^(-1) )   
  
}
