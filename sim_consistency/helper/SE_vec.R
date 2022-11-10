
SE_vec = function(j, x_eval){
  return( (shapley(j, x_eval) - shapley_popul(j, x_eval))^2 )
}


