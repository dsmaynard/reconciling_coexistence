


rm(list=ls())

set.seed(10)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

source("functions_intertidal.R")
source("functions_algae.R")

nspp <- 10

P <- matrix(runif(nspp^2),nspp,nspp)

P <- P-t(P)


x_equil <- rexp(nspp)
x_equil <- x_equil/sum(x_equil)

result_zero_sum <- find_closest_matrix_weighted(P=P, x=x_equil)


Pfit <- convert_sol_to_matrix(result_zero_sum$X,nspp)


plot(Pfit~P)

# check to make sure it's a solution
Pfit%*%x_equil

# create a weight matrix, for example, reflecting different sample sizes for each entry
weight_mat <- matrix(runif(nspp^2),nspp,nspp)

result_zero_sum_w <- find_closest_matrix_weighted(P=P, x=x_equil, weight_mat = weight_mat)


Pfit_w <- convert_sol_to_matrix(result_zero_sum_w$X,nspp)


plot(Pfit~Pfit_w)

plot(P~Pfit_w)

# check to make sure it's a solution
Pfit_w%*%x_equil



# generate some growth rates
r <- runif(nspp)


# generate an interaction matrix
A <- -matrix(runif(nspp^2), nspp,nspp)
diag(A) <- diag(A)*2

# generate an abundance vector
x_obs <- runif(nspp)

# get the best fitting result
result_LV <- fit_qp_LV(A=A,r=r,x_obs=x_obs,tol=1000)
					  
Afit <- t(matrix(result_LV$X[1:nspp^2], nspp,nspp))

rfit <- result_LV$X[(nspp^2+1):(nspp^2+nspp)]

plot(Afit~A)
plot(rfit~r)

# check to make sure it's a solution
x_obs*(rfit+Afit%*%x_obs)



# constrain the entries to be within 100% of the observed
result_LV_100 <- fit_qp_LV(A=A,r=r,x_obs=x_obs,tol=1)
					  
Afit_100 <- t(matrix(result_LV_100$X[1:nspp^2], nspp,nspp))

rfit_100 <- result_LV_100$X[(nspp^2+1):(nspp^2+nspp)]

# compare to the original. Note that it sets a bunch of interactions to be zero
plot(Afit_100~A)

# growth rate is just scales by a constant
plot(rfit_100~r)

# compare the fits
plot(Afit_100~Afit)

# check to make sure it's a solution
x_obs*(rfit_100+Afit_100%*%x_obs)

