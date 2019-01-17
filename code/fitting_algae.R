
rm(list=ls())

# !diagnostics off

# library(reshape2)
# library(reshape)
library(deSolve)
library(limSolve)
library(Matrix)
# library(viridis)
# library(RColorBrewer)
library(tidyverse)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("functions_algae.R")


# note, the units are scaled by 1/10000 relative to papers by Naughton and Narwani

A <- read_csv("../data/algae_interaction_matrix.csv") %>% as.matrix
r <- read_csv("../data/algae_growth_rates.csv") %>% unlist %>% as.numeric
x_obs <- read_csv("../data/algae_equil_abundances.csv") %>% unlist %>% as.numeric




############################# FIT THE MODEL

# fit a range of tolerance values, incremented by 10%
good_tols <- find_tolerance_qp_LV(r=r, A=A, x_obs=x_obs, tol_seq= seq(0.1,2,by=0.1))

# get the minimum tolerance value and get the predicted A
min_tol <- good_tols$tol[!good_tols$error][1]

# fit the parameters and extract
fit_pars <- fit_qp_LV(A=A,r=r,x_obs=x_obs,tol=min_tol)$X
A_fit <- t(matrix(fit_pars[1:16],4,4))
r_fit <- fit_pars[17:20]
colnames(A_fit) <- rownames(A_fit) <- names(x_obs)


# compare fits
plot(A_fit~A)
plot(r_fit~r)


#############################################
################## Plot the dynamics

# starting conditions
x0 <- rep(0.2,4) #rep(200/scale_density,nalive)


# predicted dynamics using best fit A
ts <- ode(y=x0, parms = list(r=r_fit,A=A_fit,THRESH=0), times=seq(0,500,length=500),func = LV_dyn, method="ode45")


# plot the time series
ggplot(ts %>% as.data.frame() %>% gather(species,abund,2:ncol(ts)),aes(x=time,y=abund,color=species))+geom_line()+scale_y_log10()

# compare endpoints
plot(as.numeric(ts[nrow(ts),-1])~x_obs)



# original dynamics using raw A
ts_obs <- ode(y=x0, parms = list(r=r,A=A,THRESH=0), times=seq(0,100,length=1000),func = LV_dyn, method="ode45")

# plot the time series
ggplot(ts_obs %>% as.data.frame() %>% gather(species,abund,2:ncol(ts)),aes(x=time,y=abund,color=species))+geom_line()+scale_y_log10()

# compare endpoints
plot(as.numeric(ts_obs[nrow(ts),-1])~x_obs)



