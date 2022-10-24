# start function: abh. von X,Y, subsets

model_subset = function(X, Y, subset_names){
  dat = data.frame(Y,X)
  Xnam = subset_names
 


  
  model.npp = npreg(reformulate(Xnam, "Y"), regtype = "ll", data = dat, 
                    cfac.init=1, 
                    initc.dir=1, nmulti=0,
                    remin=FALSE, bwscaling=FALSE, bwmethod="cv.aic")
  

  return(model.npp)
}



