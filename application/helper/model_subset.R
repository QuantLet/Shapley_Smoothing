# start function: abh. von X,Y, subsets

model_subset = function(X, Y, subset_names, alt, bw){
  dat = data.frame(Y,X)
  Xnam = subset_names
 

  if(alt==FALSE){
  model.final = npreg(reformulate(Xnam, "Y"), regtype = "ll", data = dat, bws=bw)
    
    
  } else{
  model.final = randomForest(reformulate(Xnam, "Y"), data = dat, ntree=500)
  }
  


  
  return(model.final)
}



