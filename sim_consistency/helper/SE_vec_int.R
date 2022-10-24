
#scalar:
#source("shapley.R")
#source("shapley_popul.R")


SE_vec_int = function(j, x_eval){
  return( (shapley_int(j, x_eval) - shapley_popul(j, x_eval))^2 ) #gerade mit skalar input
}


