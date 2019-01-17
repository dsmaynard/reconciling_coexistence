rm(list=ls())

library(limSolve)
library(Rglpk)
library(deSolve)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

source("functions_intertidal.R")


# get the latin binomial names
nmatch<-rbind(data.frame(read_csv("../data/intertidal_species_names.csv")),c("Other","Other"))

# get the counts of displacements
D<-read_csv("../data/intertidal_displacement_matrix.csv") %>% data.frame()
rownames(D) <- colnames(D)


# species abundances
x_obs<-c(0.664,0.018,0.040,0.005,0.017,0.059,0.053,0.017,0.018,0.006,0.053,0.032,0.018)
names(x_obs)<-c("B","BG","CV","HAL","MT","PP","SC","fila","flesha","crusts","R","ephem","Other")
x_obs<-sort(x_obs,decreasing = T)

# get the names for later
ord_names<-names(x_obs)


# spec<-names(x_obs)
D <- as.matrix(D[ord_names,ord_names])


# get the original H and P
H <- as.matrix(D/(D+t(D)))
P <- as.matrix(H-t(H))

# calculate the equilibrium 
find_optimal_strategy(H)

# get bet solution and optimal strategies and plot
Pfit<-find_optimal_weights(D,x_obs)$P[names(x_obs),names(x_obs)] %>% as.matrix()

# fit the mussel removal best fit
Pfit_mus <- find_optimal_weights(D[-1,-1],x_obs[-1])$P[names(x_obs)[-1],names(x_obs)[-1]]
Pfit_mus <- cbind(0,rbind(0,Pfit_mus))
colnames(Pfit_mus) <- rownames(Pfit_mus) <- names(x_obs)

### compare

plot(as.matrix(Pfit)~as.matrix(P))



# get the time series for the observed data
obs_dyn <- read_csv("../data/time_series_all.csv") %>% select(time,ord_names) %>% as.matrix 



# get the initial abundance at time=0
start<-as.numeric(obs_dyn[1,-1])
start<-as.numeric(start/sum(start))


# use the best fitting scaling
dyn_emp<-integrate_tidal_dynamics(x0=start, pars=list(P=P), maxtime=100, nsteps=34)
	
dyn_fit<-integrate_tidal_dynamics(x0=start,	pars=list(P=Pfit), maxtime=212, nsteps=34)


ggplot(dyn_emp %>% data.frame() %>%  gather("species","abundance",-time),aes(x=time,y=abundance, color=species)) + geom_line()+theme_bw()

ggplot(dyn_fit %>% data.frame() %>%  gather("species","abundance",-time),aes(x=time,y=abundance, color=species)) + geom_line()+theme_bw()


# read in the time series
obs_muss_dyn <- read_csv("../data/time_series_muss_removal.csv") %>% select(time,ord_names) %>% as.matrix


# get the starting conditions and set the mussel's abundance to zero
start_mus<-obs_muss_dyn[1,-1]
start_mus[names(start_mus)%in%c("B")]<-0
start_mus<-start_mus/sum(start_mus)

# the o
dyn_mus_orig<-integrate_tidal_dynamics(x0=start_mus, pars=list(P=Pfit), maxtime=11, nsteps=25)
dyn_mus_remove<-integrate_tidal_dynamics(x0=start_mus, pars=list(P=Pfit_mus), maxtime=9.5, nsteps=25)

ggplot(dyn_mus_orig %>% data.frame() %>%  gather("species","abundance",-time),aes(x=time,y=abundance, color=species)) + geom_line()+theme_bw()

ggplot(dyn_mus_remove %>% data.frame() %>%  gather("species","abundance",-time),aes(x=time,y=abundance, color=species)) + geom_line()+theme_bw()


