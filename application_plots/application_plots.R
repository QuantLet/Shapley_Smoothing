setwd("/Users/ratmir/MISE/application")
source("functions.R")

cov=0
sigma_sim=matrix(c(1, cov,
                   cov, 1), nrow=2, ncol=2)

source("integral_estimation_2d.R") 
source("shapley_int.R")

library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)
#load data
f <- file.choose()

data=read.csv(file=f, header = TRUE, sep = ",", dec = ".")                  

data_imp = data[, c("price.baseMSRP","BLUETOOTH","HEATED.SEATS", "NAVIGATION", 
"features.Measurements.Curb.weight.lbs", "features.Engine.Horsepower.hp", "year", 
"features.Measurements.Width.in.", "features.Measurements.Length.in.", 
"features.Measurements.Height.in.", "Total.Seating", "features.Tires.and.Wheels.Tire.Aspect.Ratio", 
"features.Measurements.Wheel.base.in.", "features.Tires.and.Wheels.Wheel.Diameter")]
data_imp[,1] = data_imp[,1]/1000

#binary
data_imp$BLUETOOTH = ifelse(data_imp$BLUETOOTH == 'Yes', 1, 0)
data_imp$BLUETOOTH = as.factor(data_imp$BLUETOOTH)
data_imp$HEATED.SEATS = ifelse(data_imp$HEATED.SEATS == 'Yes', 1, 0)
data_imp$HEATED.SEATS = as.factor(data_imp$HEATED.SEATS)
data_imp$NAVIGATION = ifelse(data_imp$NAVIGATION == 'Yes', 1, 0)
data_imp$NAVIGATION = as.factor(data_imp$NAVIGATION)

names(data_imp) = c("price", "BLUETOOTH","heated","navi","weight", "hp", "year",
                    "width", "length", "height", "seats", "ratio", "wheels", "diameter")


#stratify
dat = data_imp[,c("price", "hp", "year", "weight", "length")]
summary(dat)

dat_sort = dat[order(dat[,"year"]),]

dat_first = dat_sort[1:12230, ]
dat_second = dat_sort[12231:24185, ]
dat_third = dat_sort[24186:38435, ]

#All possible subsets, specify time interval
X = within(dat_first, rm("year", "price"))
d = ncol(X)
names = names(X)
names(X) = c("X1", "X2", "X3") 
Y = dat_first$price

#summary statistics
par(mar = c(1.3, 2.5, 3, 0.4))
boxplot(X[!(X[,1]>400 | X[,1]<120),1], cex.main=3,
        cex.axis=2, font.main = 2, main="Horsepower")
boxplot(X[!(X[,2]<2000 | X[,2]>6000),2], cex.main=3,
        cex.axis=2, font.main = 2, main="Weight")
boxplot(X[!(X[,3]<170 | X[,3]>260),3], cex.main=3,
        cex.axis=2, font.main = 2, main="Length")

subs <<- subsets(X)

#Get model fits and sort them in a list
sub_bw = subs

c=3  
sub_bw[[1]][,1] = 1*c
sub_bw[[1]][,2] = 40*c
sub_bw[[1]][,3] = 2*c # ok

sub_bw[[2]][1,1] = 4*c
sub_bw[[2]][2,1] = 60*c

sub_bw[[2]][1,2] = 24*c
sub_bw[[2]][2,2] = 6*c

sub_bw[[2]][1,3] = 80*c
sub_bw[[2]][2,3] = 1*c # ?

sub_bw[[3]][1,1] = 12*c
sub_bw[[3]][2,1] = 80*c
sub_bw[[3]][3,1] = 4*c

model_list <<- model_list_fct(subs=subs, alt=FALSE, sub_bw) 


grid_a = seq(120, 400, length.out=100)
grid_b = seq(2000, 6000, length.out=100)
grid_c = seq(170, 260, length.out=100)

  
eval_a = rbind(grid_a, rep(3500, 100),  rep(190, 100))
col_a = shapley_vec(j=1, x_eval=eval_a , alt=FALSE, model_list = model_list)

eval_b = rbind(rep(250, 100), grid_b, rep(190, 100))
col_b = shapley_vec(j=2, x_eval=eval_b , alt=FALSE, model_list = model_list)

eval_c = rbind(rep(250, 100), rep(3500, 100), grid_c) 
col_c = shapley_vec(j=3, x_eval=eval_c, alt=FALSE, model_list = model_list)

#Sequential plots over x1:
eval_a = rbind(grid_a, rep(3500, 100),  rep(190, 100))
col_a = shapley_vec(j=1, x_eval=eval_a , alt=FALSE, model_list = model_list)
col_b = shapley_vec(j=2, x_eval=eval_a , alt=FALSE, model_list = model_list)
col_c = shapley_vec(j=3, x_eval=eval_a, alt=FALSE, model_list = model_list)

#Sequential plots over x2:
eval_b = rbind(rep(250, 100), grid_b,  rep(190, 100))
col_a = shapley_vec(j=1, x_eval=eval_b , alt=FALSE, model_list = model_list)
col_b = shapley_vec(j=2, x_eval=eval_b , alt=FALSE, model_list = model_list)
col_c = shapley_vec(j=3, x_eval=eval_b, alt=FALSE, model_list = model_list)

#Sequential plots over x3:
eval_c = rbind(rep(250, 100), rep(3500, 100),  grid_c)
col_a = shapley_vec(j=1, x_eval=eval_c , alt=FALSE, model_list = model_list)
col_b = shapley_vec(j=2, x_eval=eval_c , alt=FALSE, model_list = model_list)
col_c = shapley_vec(j=3, x_eval=eval_c, alt=FALSE, model_list = model_list)


#Get prediction
pred_fct = function(eval){
x_eval = t(eval)
x_eval_vec <- data.frame(X1 = c(x_eval[,1]), X2 = c(x_eval[,2]),  X3 = c(x_eval[,3]))  
names_vars = model_list[[7]]$xnames 
pred = predict(model_list[[7]], newdata = x_eval_vec[, c(names_vars), drop = FALSE ]) 
return(pred)
}

pred = pred_fct(eval = eval_c)


#Slice plots with CI
par(mar = c(4.2, 2.2, 0.7, 0.4)) # for overleaf 620x498 
plot(y=col_a, x=grid_a, type="l", xlab="Horsepower [hp], 2014 - 2020", ylab="", ylim=c(-20,70), main="", cex.main=1.5,
     cex.axis=1.5, font.main = 1.5, cex.lab=1.5)

x1_ci = seq(120, 400, length.out=48)
points(x=x1_ci, y=points[1,] ,cex=0.3, col="red")
points(x=x1_ci, y=points[2,] ,cex=0.3, col="red")
lines(x=x1_ci, y=points[1,] ,cex=0.3, col="red")
lines(x=x1_ci, y=points[2,] ,cex=0.3, col="red")
abline(a=0,b=0, lty=2)


plot(density(X[,1]), xlab="Horsepower [hp]", ylab="",main="", cex.main=1.5,
     cex.axis=1.5, font.main = 1.5, cex.lab=1.5)

plot(density(X[,2]), xlab="Weight [lbs]", ylab="",main="", cex.main=1.5,
     cex.axis=1.5, font.main = 1.5, cex.lab=1.5)

plot(density(X[,3]), xlab="Length [in]", ylab="",main="", cex.main=1.5,
     cex.axis=1.5, font.main = 1.5, cex.lab=1.5)


points(x=x1_ci, y=points[7,] ,cex=0.3, col="blue")
points(x=x1_ci, y=points[8,] ,cex=0.3, col="blue")
lines(x=x1_ci, y=points[7,] ,cex=0.3, col="blue")
lines(x=x1_ci, y=points[8,] ,cex=0.3, col="blue")

plot(y=col_b, x=grid_b, type="l", xlab="Weight [lbs], 2014 - 2020", ylab="", ylim=c(-35,30), main="", cex.main=1.5,
     cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
x2_ci = seq(2000, 6000, length.out=48)
points(x=x2_ci, y=points[3,] ,cex=0.3, col="red")
points(x=x2_ci, y=points[4,] ,cex=0.3, col="red")
lines(x=x2_ci, y=points[3,] ,cex=0.3, col="red")
lines(x=x2_ci, y=points[4,] ,cex=0.3, col="red")
abline(a=0,b=0, lty=2)


points(x=x2_ci, y=points[9,] ,cex=0.3, col="blue")
points(x=x2_ci, y=points[10,] ,cex=0.3, col="blue")
lines(x=x2_ci, y=points[9,] ,cex=0.3, col="blue")
lines(x=x2_ci, y=points[10,] ,cex=0.3, col="blue")

plot(y=col_c, x=grid_c, type="l", xlab="Length [in], 2014 - 2020", ylab="", ylim=c(-35,14), main="", cex.main=1.5,
     cex.axis=1.5, font.main = 1.5, cex.lab=1.5)
x3_ci = seq(170, 260, length.out=48)
points(x=x3_ci, y=points[5,] ,cex=0.3, col="red")
points(x=x3_ci, y=points[6,] ,cex=0.3, col="red")
lines(x=x3_ci, y=points[5,] ,cex=0.3, col="red")
lines(x=x3_ci, y=points[6,] ,cex=0.3, col="red")
abline(a=0,b=0, lty=2)


points(x=x3_ci, y=points[11,] ,cex=0.3, col="blue")
points(x=x3_ci, y=points[12,] ,cex=0.3, col="blue")
lines(x=x3_ci, y=points[11,] ,cex=0.3, col="blue")
lines(x=x3_ci, y=points[12,] ,cex=0.3, col="blue")


#3D plots

x1_grid = seq(150, 350, length.out=100) 
x2_grid = seq(2000, 6000, length.out=100)


grid = t(expand.grid(x1_grid, x2_grid))
grid = rbind(grid , rep(190,ncol(grid)) )

shap_eval_est1=shapley_vec(j=1, grid, alt=FALSE, model_list = model_list) 
shap_eval_est2=shapley_vec(j=2, grid, alt=FALSE, model_list = model_list)

#On data points
#shap_eval_est1=shapley_vec(j=1, t(X), alt=TRUE, model_list = model_list) 
#shap_eval_est2=shapley_vec(j=2, t(X), alt=TRUE, model_list = model_list)

#plot(y=shap_eval_est1, x=X[,1])
#plot(y=shap_eval_est2, x=X[,2])



#exclude outliers

outl_one = shap_eval_est1<=100 & shap_eval_est1>=-100
outl_two = shap_eval_est2<=100 & shap_eval_est2>=-100


shap_eval_est1[!outl_one] = median(shap_eval_est1)
shap_eval_est2[!outl_two] = median(shap_eval_est2)





#On X3:

x2_grid = seq(2000, 6000, length.out=100) 
x3_grid = seq(170, 260, length.out=100)

#grid = t(expand.grid(x1_grid, x3_grid))
#grid = rbind(rep(250,ncol(grid)),grid)

grid = t(expand.grid(x2_grid, x3_grid))
grid = rbind(rep(250,ncol(grid)), grid[1,], grid[2,])


shap_eval_est3=shapley_vec(j=3, grid, alt=FALSE, model_list = model_list) 


surface3_est=t(pracma::Reshape(shap_eval_est3, length(x2_grid), length(x3_grid))) 





surface1_est=t(pracma::Reshape(shap_eval_est1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2_est=t(pracma::Reshape(shap_eval_est2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 



library(plotly)
par(mar = c(0, 0, 0, 0))
fig1_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_est, type="surface") %>% hide_colorbar()
fig1_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig2_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_est, type="surface") %>% hide_colorbar()
fig2_est <- fig2_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig3_est = plot_ly(x=sort(x2_grid), y=sort(x3_grid), z=surface3_est, type="surface") %>% hide_colorbar()
fig3_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x3'),zaxis=list(title='shap3')))


mrg <- list(l = 0, r = 0,
            b = 0, t = 0,
            pad = 0)

#different colors
plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface1_est, opacity = 0.8) %>%
  layout(scene=list(xaxis=list(title='Horsepower [hp]'), yaxis=list(title='Weight [lbs]'), zaxis=list(title='', range = c(-25,70)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)

#different colors
plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface2_est, opacity = 0.8) %>%
  layout(scene=list(xaxis=list(title='Horsepower [hp]'), yaxis=list(title='Weight [lbs]'), zaxis=list(title='', range = c(-15,40)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)

#different colors
plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x2_grid), y=~sort(x3_grid),z = ~surface3_est, opacity = 0.8) %>%
  layout(scene=list(xaxis=list(title='Horsepower [hp]'), yaxis=list(title='Length [in]'), zaxis=list(title='', range = c(-20,8)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)




