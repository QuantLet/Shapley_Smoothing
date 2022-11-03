model_list_fct = function(subs, alt, sub_bw){
names = colnames(X)
model_list_outer = vector(mode="list")


for (j in 1:ncol(X)){
  model_list_inner = vector(mode="list")
  for (i in 1:ncol(subs[[j]]) ){
    
    sub_nam = names[subs[[j]][,i]]
    bw = sub_bw[[j]][,i]
    model_list_inner[[i]] = model_subset(X=X, Y=Y, subset_names = sub_nam, alt=alt, bw) 
  }
  model_list_outer[[j]] = model_list_inner
}
#list containing model fits
model_l = unlist(model_list_outer, recursive = FALSE)


return(model_l)
}