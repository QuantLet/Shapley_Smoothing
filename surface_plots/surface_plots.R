#Load DGP

library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)

source("functions.R") # load all functions

#define cov before reading the functions
cov=0
sigma_sim=matrix(c(4, cov, cov,
                   cov, 4, cov,
                   cov, cov, 4), nrow=3, ncol=3)


source("integral_population.R") #adjust DGP here also.
source("integral_estimation.R")
source("shapley_int.R")
source("SE_vec_int.R")


g1 = function(X){ return( -sin(2*X[,1]) ) } 
g2 = function(X){ return( cos(3*X[,2])  ) }  
g3 = function(X){ return( 0.5*X[,3] ) } 
int = function(X){
  x1 = X[,1]
  x2 = X[,2]
  return( 2*cos(x1)*sin(2*x2)  ) 
}

m_full_why = function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return(-sin(2*x1) + cos(3*x2) + 0.5*x3 + 2*cos(x1)*sin(2*x2)  )
}


true_model_list = list()
true_model_list[[1]] = m_x1
true_model_list[[2]] = m_x2
true_model_list[[3]] = m_x3 
true_model_list[[4]] = m_x1_x2
true_model_list[[5]] = m_x1_x3
true_model_list[[6]] = m_x2_x3 
true_model_list[[7]] = m_full_why 


l = -2; u = 2; N=1000
l_int = l; u_int = u
d = 3

X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))

#DGP
Y <<- g1(X) + g2(X) + g3(X) + int(X) + rnorm(N)


#All possible subsets
subs <<- subsets(X)

#Get model fits and sort them in a list
model_list <<- model_list_fct(subs=subs)

x1_grid = seq(-2, 2, length.out=30) 
x2_grid = seq(-2, 2, length.out=30)

grid=t(expand.grid(x1_grid, x2_grid))
grid = rbind(grid , rep(0,ncol(grid)))

#IB Shapley Curves
shap_eval_est1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval_est2 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
  shap_eval_est1[i] = shapley_int(j=1, grid_col) 
  shap_eval_est2[i] = shapley_int(j=2, grid_col) 
  
}


#CB shapley curves
shap_eval_est1=shapley_vec(j=1, grid) 
shap_eval_est2=shapley_vec(j=2, grid)


#population shapley curves
shap_eval1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval2 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
  shap_eval1[i] = shapley_popul(j=1, grid_col) 
  shap_eval2[i] = shapley_popul(j=2, grid_col) 
  
}



shap_SE1 = SE_vec(j=1, grid)
surface1_SE=t(pracma::Reshape(shap_SE1, length(x1_grid), length(x2_grid))) 

shap_SE2 = SE_vec(j=2, grid)
surface2_SE=t(pracma::Reshape(shap_SE2, length(x1_grid), length(x2_grid))) 


surface1=t(pracma::Reshape(shap_eval1, length(x1_grid), length(x2_grid))) 
surface2=t(pracma::Reshape(shap_eval2, length(x1_grid), length(x2_grid)))
surface1_est=t(pracma::Reshape(shap_eval_est1, length(x1_grid), length(x2_grid))) 
surface2_est=t(pracma::Reshape(shap_eval_est2, length(x1_grid), length(x2_grid))) 

surface1_SE = (surface1 - surface1_est)^2
surface2_SE = (surface2 - surface2_est)^2


library(plotly)
par(mar = c(0, 0, 0, 0))
fig1 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1, type="surface") %>% hide_colorbar()
fig1 <- fig1 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig1_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_est, type="surface") %>% hide_colorbar()
fig1_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig2 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2, type="surface") %>% hide_colorbar()
fig2 <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig2_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_est, type="surface") %>% hide_colorbar()
fig2_est <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig1_SE = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_SE, type="surface") %>% hide_colorbar()
fig1_SE <- fig1_SE %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='SE shap1')))

fig2_SE = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_SE, type="surface") %>% hide_colorbar()
fig2_SE <- fig2_SE %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='SE shap2')))


#curves plots
plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface1 , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface1_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    ))

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface2 , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface2_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))      ))



#SE plots

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface1_SE , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title='', range=c(0,1)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    )) 

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface2_SE , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title='', range=c(0,1)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    ))



