install.packages("remotes")
remotes::install_github("reyesem/IntroAnalysis")

bootstrap = function(j, B=800){

source("functions.R") 
  
subs <<- subsets(X)
  
library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)
library(IntroAnalysis)

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

#1:This is phi_j
shap_estimated1 = shapley_vec(j=1, c(ci_points[j], 3500, 190), alt=FALSE, model_list = model_list) 
shap_estimated2 = shapley_vec(j=2, c(250, ci_points2[j], 190), alt=FALSE, model_list = model_list) 
shap_estimated3 = shapley_vec(j=3, c(250 , 3500, ci_points3[j]), alt=FALSE, model_list = model_list)

#2: Calculate g as a fct of b. 
n <<- nrow(X)
g1 <<- b1*log(log(n))*0.7
g2 <<- b2*log(log(n))*0.7
g3 <<- b3*log(log(n))*0.7
g12 <<- b12*log(log(n))*0.7
g13 <<- b13*log(log(n))*0.7
g23 <<- b23*log(log(n))*0.7
g123 <<- b123*log(log(n))*0.7

model_list_bs = vector(mode="list")
model.np1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = X, bws=b1)
model.np2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = X, bws=b2) 
model.np3 = npreg(reformulate("X3", "Y"), regtype = "ll", data = X, bws=b3) 

model.np12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = X, bws=b12) 
model.np13 = npreg(reformulate(c("X1", "X3"), "Y"), regtype = "ll", data = X, bws=b13) 
model.np23 = npreg(reformulate(c("X2", "X3"), "Y"), regtype = "ll", data = X, bws=b23) 

model.np123 = npreg(reformulate(c("X1", "X2", "X3"), "Y"), regtype = "ll", data = X, bws=b123) 

model_list_bs[[1]] = model.np1 
model_list_bs[[2]] = model.np2 
model_list_bs[[3]] = model.np3 
model_list_bs[[4]] = model.np12 
model_list_bs[[5]] = model.np13 
model_list_bs[[6]] = model.np23 
model_list_bs[[7]] = model.np123 

#3:calculate the plain residuals based on bandwidth b
eps1 = Y - predict(model_list_bs[[1]], newdata = X) 
eps2 = Y - predict(model_list_bs[[2]], newdata = X) 
eps3 = Y - predict(model_list_bs[[3]], newdata = X) 
eps12 = Y - predict(model_list_bs[[4]], newdata = X) 
eps13 = Y - predict(model_list_bs[[5]], newdata = X) 
eps23 = Y - predict(model_list_bs[[6]], newdata = X) 
eps123 = Y - predict(model_list_bs[[7]], newdata = X) 

df1_pred = data.frame(X1 = c(ci_points[j]), X2 = c(3500), X3=c(190) )
df2_pred = data.frame(X1 = c(250), X2 = c(ci_points2[j]), X3=c(190) )
df3_pred = data.frame(X1 = c(250), X2 = c(3500), X3=c(ci_points3[j]) )

mod_g1 = npreg(reformulate(c("X1"), "Y"), regtype = "ll", data = X, bws=g1) 
mod_g2 = npreg(reformulate(c("X2"), "Y"), regtype = "ll", data = X, bws=g2) 
mod_g3 = npreg(reformulate(c("X3"), "Y"), regtype = "ll", data = X, bws=g3) 
mod_g12 = npreg(reformulate(c("X1","X2"), "Y"), regtype = "ll", data = X, bws=g12) 
mod_g13 = npreg(reformulate(c("X1","X3"), "Y"), regtype = "ll", data = X, bws=g13) 
mod_g23 = npreg(reformulate(c("X2", "X3"), "Y"), regtype = "ll", data = X, bws=g23) 
mod_g123 = npreg(reformulate(c("X1", "X2", "X3"), "Y"), regtype = "ll", data = X, bws=g123) 

pred1 = predict(mod_g1, newdata = df1_pred)
pred2 = predict(mod_g2, newdata = df1_pred)
pred3 = predict(mod_g3, newdata = df1_pred)
pred12 = predict(mod_g12, newdata = df1_pred)
pred13 = predict(mod_g13, newdata = df1_pred)
pred23 = predict(mod_g23, newdata = df1_pred)
pred123 = predict(mod_g123, newdata = df1_pred)

pred1_2 = predict(mod_g1, newdata = df2_pred)
pred2_2 = predict(mod_g2, newdata = df2_pred)
pred3_2 = predict(mod_g3, newdata = df2_pred)
pred12_2 = predict(mod_g12, newdata = df2_pred)
pred13_2 = predict(mod_g13, newdata = df2_pred)
pred23_2 = predict(mod_g23, newdata = df2_pred)
pred123_2 = predict(mod_g123, newdata = df2_pred)

pred1_3 = predict(mod_g1, newdata = df3_pred)
pred2_3 = predict(mod_g2, newdata = df3_pred)
pred3_3 = predict(mod_g3, newdata = df3_pred)
pred12_3 = predict(mod_g12, newdata = df3_pred)
pred13_3 = predict(mod_g13, newdata = df3_pred)
pred23_3 = predict(mod_g23, newdata = df3_pred)
pred123_3 = predict(mod_g123, newdata = df3_pred)
                      
for ( i in 1:B){

  
  #4: calculate wild residuals, just multiply with V
  random = rmammen(n=n, construct= c("two-point mass"))
  
  #5
  Y_g1 = pred1 + eps1*random
  Y_g2 = pred2 + eps2*random  
  Y_g3 = pred3 + eps3*random  
  Y_g12 = pred12 + eps12*random  
  Y_g13 = pred13 + eps13*random  
  Y_g23 = pred23 + eps23*random  
  Y_g123 = pred123 + eps123*random  
  
  #6: predict with bandwidth b on (new Y, old X) 
  model.1 = npreg(reformulate(c("X1"), "Y_g1"), regtype = "ll", data = data.frame(Y_g1, X), bws=b1) 
  model.2 = npreg(reformulate(c("X2"), "Y_g2"), regtype = "ll", data = data.frame(Y_g2, X), bws=b2) 
  model.3 = npreg(reformulate(c("X3"), "Y_g3"), regtype = "ll", data = data.frame(Y_g3, X), bws=b3) 
  model.12 = npreg(reformulate(c("X1", "X2"), "Y_g12"), regtype = "ll", data = data.frame(Y_g12, X), bws=b12) 
  model.13 = npreg(reformulate(c("X1", "X3"), "Y_g13"), regtype = "ll", data = data.frame(Y_g13, X), bws=b13) 
  model.23 = npreg(reformulate(c("X2", "X3"), "Y_g23"), regtype = "ll", data = data.frame(Y_g23, X), bws=b23) 
  model.123 = npreg(reformulate(c("X1", "X2", "X3"), "Y_g123"), regtype = "ll", data = data.frame(Y_g123, X), bws=b123) 
  
  pred_bs1 = predict(model.1, newdata = df1_pred)
  pred_bs2 = predict(model.2, newdata = df1_pred)
  pred_bs3 = predict(model.3, newdata = df1_pred)
  pred_bs12 = predict(model.12, newdata = df1_pred)
  pred_bs13 = predict(model.13, newdata = df1_pred)
  pred_bs23 = predict(model.23, newdata = df1_pred)
  pred_bs123 = predict(model.123, newdata = df1_pred)
  
  pred_bs1_2 = predict(model.1, newdata = df2_pred)
  pred_bs2_2 = predict(model.2, newdata = df2_pred)
  pred_bs3_2 = predict(model.3, newdata = df2_pred)
  pred_bs12_2 = predict(model.12, newdata = df2_pred)
  pred_bs13_2 = predict(model.13, newdata = df2_pred)
  pred_bs23_2 = predict(model.23, newdata = df2_pred)
  pred_bs123_2 = predict(model.123, newdata = df2_pred)
  
  pred_bs1_3 = predict(model.1, newdata = df3_pred)
  pred_bs2_3 = predict(model.2, newdata = df3_pred)
  pred_bs3_3 = predict(model.3, newdata = df3_pred)
  pred_bs12_3 = predict(model.12, newdata = df3_pred)
  pred_bs13_3 = predict(model.13, newdata = df3_pred)
  pred_bs23_3 = predict(model.23, newdata = df3_pred)
  pred_bs123_3 = predict(model.123, newdata = df3_pred)

  #7:Final step: BS SHAP
  bs_shap[i,1] = 0.33*(pred_bs1 - pred1) - 0.166*(pred_bs2 - pred2) - 0.166*(pred_bs3 - pred3) + 0.166*(pred_bs12 - pred12) + 0.166*(pred_bs13 - pred13)-
   0.33*(pred_bs23 - pred23) + 0.33*(pred_bs123 - pred123) -0.3333*(mean(Y_g123) - mean(Y))

  bs_shap[i,2] = -0.166*(pred_bs1_2 - pred1_2) + 0.33*(pred_bs2_2 - pred2_2) - 0.166*(pred_bs3_2 - pred3_2) + 0.166*(pred_bs12_2 - pred12_2) + 0.166*(pred_bs13_2 - pred13_2)-
    0.33*(pred_bs23_2 - pred23_2) + 0.33*(pred_bs123_2 - pred123_2) -0.3333*(mean(Y_g123) - mean(Y))
  
  bs_shap[i,3] = -0.166*(pred_bs1_3 - pred1_3) -0.166*(pred_bs2_3 - pred2_3) + 0.33*(pred_bs3_3 - pred3_3) - 0.33*(pred_bs12_3 - pred12_3) + 0.166*(pred_bs13_3 - pred13_3)+
    0.166*(pred_bs23_3 - pred23_3) + 0.33*(pred_bs123_3 - pred123_3) -0.3333*(mean(Y_g123) - mean(Y))
  
  
  
}

alpha = 0.05

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


return(c(qs1, qs2, qs3))
}


source("functions.R") 

#load data
f <- file.choose() 

data=read.csv(file=f,               
              header = TRUE,        
              sep = ",",            
              dec = ".")                  

data_imp = data[, c("price.baseMSRP","BLUETOOTH","HEATED.SEATS", "NAVIGATION", 
                    "features.Measurements.Curb.weight.lbs", "features.Engine.Horsepower.hp", "year", 
                    "features.Measurements.Width.in.", "features.Measurements.Length.in.", 
                    "features.Measurements.Height.in.", "Total.Seating", "features.Tires.and.Wheels.Tire.Aspect.Ratio", 
                    "features.Measurements.Wheel.base.in.", "features.Tires.and.Wheels.Wheel.Diameter")]

data_imp[,1] = data_imp[,1]/1000
names(data_imp) = c("price", "BLUETOOTH","heated","navi","weight", "hp", "year",
                    "width", "length", "height", "seats", "ratio","wheels", "diameter")

dat = data_imp[,c("price", "hp", "year", "weight", "length")]

dat_sort = dat[order(dat[,"year"]),]
dat_first = dat_sort[1:12230, ]  
#dat_second = dat_sort[12231:24185, ]
#dat_third = dat_sort[24186:38435, ]

#All possible subsets
X <<- within(dat_first, rm("year", "price"))
d <<- ncol(X)
names = names(X)
names(X) = c("X1", "X2", "X3") 
Y <<- dat_first$price

dat_bs <<- cbind(Y,X)


library(parallel)
qs = mclapply(1:48, bootstrap, mc.cores=48) # 48 grid points



points=matrix(unlist(qs), byrow = FALSE, ncol=48)
# 1+2 ->first variable; 3+4 -> second variable; 5+6 -> third variable



