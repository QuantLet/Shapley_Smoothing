
bootstrap = function(j, B=1000){

source("functions.R") 
  
  
subs <<- subsets(X)
  
  

library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)
library(parallel)

sub_bw = subs; c=3
#d=3
sub_bw[[1]][,1] = 1*c
sub_bw[[1]][,2] = 40*c
sub_bw[[1]][,3] = 2*c

sub_bw[[2]][1,1] = 4*c
sub_bw[[2]][2,1] = 60*c

sub_bw[[2]][1,2] = 24*c
sub_bw[[2]][2,2] = 6*c

sub_bw[[2]][1,3] = 80*c
sub_bw[[2]][2,3] = 1*c

sub_bw[[3]][1,1] = 12*c
sub_bw[[3]][2,1] = 80*c
sub_bw[[3]][3,1] = 4*c

model_list <<- model_list_fct(subs=subs, alt=FALSE, sub_bw) 



#generate B bootstrap samples
bs_shap = bs_shap=matrix(0, ncol=3, nrow=B)

#take bandwidths
b1 <<- model_list[[1]]$bw
b2 <<- model_list[[2]]$bw
b3 <<- model_list[[3]]$bw
b12 <<- model_list[[4]]$bw
b13 <<- model_list[[5]]$bw
b23 <<- model_list[[6]]$bw
b123 <<- model_list[[7]]$bw



ci_points = seq(120, 400, length.out=48)
ci_points2 = seq(2000, 6000, length.out=48)
ci_points3 = seq(170, 260, length.out=48)


shap_estimated1 = shapley_vec(j=1, c(ci_points[j], 3500, 190), alt=FALSE, model_list = model_list) 
shap_estimated2 = shapley_vec(j=2, c(250, ci_points2[j], 190), alt=FALSE, model_list = model_list) 
shap_estimated3 = shapley_vec(j=3, c(250 , 3500, ci_points3[j]), alt=FALSE, model_list = model_list)

for ( i in 1:B){
  model_list_bs = vector(mode="list")
  samp = dat_bs[sample(nrow(dat_bs), 11955, replace=TRUE), ]
  model.np1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = samp, bws=b1)
  model.np2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = samp, bws=b2) 
  model.np3 = npreg(reformulate("X3", "Y"), regtype = "ll", data = samp, bws=b3) 
  
  model.np12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = samp, bws=b12) 
  model.np13 = npreg(reformulate(c("X1", "X3"), "Y"), regtype = "ll", data = samp, bws=b13) 
  model.np23 = npreg(reformulate(c("X2", "X3"), "Y"), regtype = "ll", data = samp, bws=b23) 
  
  model.np123 = npreg(reformulate(c("X1", "X2", "X3"), "Y"), regtype = "ll", data = samp, bws=b123) 
  
  
  model_list_bs[[1]] = model.np1 
  model_list_bs[[2]] = model.np2 
  model_list_bs[[3]] = model.np3 
  model_list_bs[[4]] = model.np12 
  model_list_bs[[5]] = model.np13 
  model_list_bs[[6]] = model.np23 
  model_list_bs[[7]] = model.np123 
  
  
  bs_shap[i,1] = shapley_vec(j=1, c(ci_points[j], 3500, 190), alt=FALSE, model_list = model_list_bs) 
  bs_shap[i,2] = shapley_vec(j=2, c(250, ci_points2[j], 190), alt=FALSE, model_list = model_list_bs) 
  bs_shap[i,3] = shapley_vec(j=3, c(250 , 3500, ci_points3[j]), alt=FALSE, model_list = model_list_bs) 
  
}

alpha = 0.05
qs1q = quantile(bs_shap[,1], probs=c(alpha/2, 1 - alpha/2) )
qs2q = quantile(bs_shap[,2], probs=c(alpha/2, 1 - alpha/2) )
qs3q = quantile(bs_shap[,3], probs=c(alpha/2, 1 - alpha/2) )

# SE approach
sd1 = sd(bs_shap[,1])
CI_lower1 = shap_estimated1 - qnorm(1 - alpha / 2)*sd1
CI_upper1 = shap_estimated1 + qnorm(1 - alpha / 2)*sd1
qs1 = c(CI_lower1, CI_upper1)

sd2 = sd(bs_shap[,2])
CI_lower2 = shap_estimated2 - qnorm(1 - alpha / 2)*sd2
CI_upper2 = shap_estimated2 + qnorm(1 - alpha / 2)*sd2
qs2 = c(CI_lower2, CI_upper2)

sd3 = sd(bs_shap[,3])
CI_lower3 = shap_estimated3 - qnorm(1 - alpha / 2)*sd3
CI_upper3 = shap_estimated3 + qnorm(1 - alpha / 2)*sd3
qs3 = c(CI_lower3, CI_upper3)




return(c(qs1, qs2, qs3, qs1q, qs2q, qs3q))
}


source("functions.R") 

f <- file.choose() 

data=read.csv(file=f, header = TRUE, sep = ",", dec = ".")                  

data_imp = data[, c("price.baseMSRP","BLUETOOTH","HEATED.SEATS", "NAVIGATION", 
                    "features.Measurements.Curb.weight.lbs", "features.Engine.Horsepower.hp", "year", 
                    "features.Measurements.Width.in.", "features.Measurements.Length.in.", 
                    "features.Measurements.Height.in.", "Total.Seating", "features.Tires.and.Wheels.Tire.Aspect.Ratio", 
                    "features.Measurements.Wheel.base.in.", "features.Tires.and.Wheels.Wheel.Diameter")]

data_imp[,1] = data_imp[,1]/1000
names(data_imp) = c("price", "BLUETOOTH","heated","navi","weight", "hp", "year",
                    "width", "length", "height", "seats", "ratio","wheels", "diameter")

#stratify
dat = data_imp[,c("price", "hp", "year", "weight", "length")]

dat_sort = dat[order(dat[,"year"]),]
#dat_first = dat_sort[1:12230, ]  
dat_second = dat_sort[12231:24185, ]
#dat_third = dat_sort[24186:38435, ]

X <<- within(dat_second, rm("year", "price"))
d <<- ncol(X)
names = names(X)
names(X) = c("X1", "X2", "X3") 
Y <<- dat_second$price
dat_bs <<- cbind(Y,X)


qs = mclapply(1:48, bootstrap, mc.cores=48)
points=matrix(unlist(qs), byrow = FALSE, ncol=48)



