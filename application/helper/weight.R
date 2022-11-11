weight = function(j,k, names_vars){
  
  indicator = as.integer(sum(as.integer(paste("X",j, sep="") == c(names_vars)))>0) # j=variable, k= subset s
  sign = ifelse(indicator > 0, 1, -1) # sign fct
  
  card_s = length(names_vars)
    
  return( sign*(1/d)*( nchoosek(d-1, card_s-indicator) )^(-1) )   
  
}
