
#################################
#### Linear programming

find_optimal_strategy<- function(H,verbose=TRUE,time_limit=5000){
    n <- dim(H)[1]
    f.obj <- rep(1, n)
    In<-diag(n)
    f.con <- H
    f.rhs <- rep(1, n)
    f.dir <-rep("<=", n)
    z<-Rglpk_solve_LP(obj=f.obj,mat=f.con,dir=f.dir,rhs=f.rhs,max=TRUE,control=list(tm_limit=time_limit,presolve=TRUE,verbose=verbose))
    return(z$solution / sum(z$solution))
}


#################################
#### QUADRATIC PROGRAMMING

# given a weighting, find the best fitting zero-sum payoff matrix matching D-t(D)
find_closest_matrix_weighted <- function(D=NULL, P=NULL, x, weight_mat=NULL){
	# get H and P if not supplied
	if(!is.null(D)){
		D<-as.matrix(D)
		D<-D[names(x),names(x)]
		H<-D/(D+t(D))
		H[is.na(H)]<-0.5	
		P<-H-t(H)		
	}
	if(!is.null(weight_mat)){
		# make sure the weights are positive
		weight_mat[weight_mat==0]<-min(weight_mat[weight_mat>0],na.rm=T)/2
		# scale the weights to one
		weight_mat<-weight_mat/max(weight_mat[upper.tri(weight_mat)],na.rm=T)
	}
	x<-as.numeric(x)
	# construct the variables for quadratic programming
    n <- length(x)
    m <- choose(n, 2)
    B <- as.matrix(get_vector(P, n), ncol = 1)
    C <- get_constraints(x, n, m)
    E <- C[1:n, 1:m]
    A <- C[(n+1):(n + 2 * m), 1:m]
    Fvec <- rep(0, n)
    # creat the bounds on P. note that this is set to 0.99 rather than 1 to ensure log(H) is well defined. 
    bvec2 <- -c(rep(1, m), rep(1, m))*0.99
    # implement quadratic programming to get a solution
    if(!is.null(weight_mat)){
    	res1 <- lsei(A = diag(m), B = B, E = E, F= Fvec, G = -A, H = bvec2, Wa=get_vector(weight_mat,n))
    }
    else{
    	res1 <- lsei(A = diag(m), B = B, E = E, F= Fvec, G = -A, H = bvec2)
    }
    return (res1)
}



# search through weighted matrices and find those that give the best fit based on log-likelihood
find_optimal_weights<-function(D,x_obs){
	snames<-colnames(D)
	x_obs<-x_obs/sum(x_obs)
	n<-length(x_obs)
	# get the Bayesian variance, scaled so the max is 1
	ap<-D+1
	bp<-t(D)+1
	var_ab<-(ap*bp)/((ap+bp)^2*(ap+bp+1))
	w_list<-1/var_ab
	diag(w_list)<-NA
	w_list<-w_list/max(w_list,na.rm=T)
	# exponent for scaling the weights
	alpha<-c(1/10000000,.05,0.1,0.25,.5,1,1.5,2,3,4) 
	max_ll <- -Inf
	for(i in 1:length(alpha)){
		# raise to an exponent
		wmat<-w_list^(alpha[i])
		# make sure value isn't too small, or routine won't converge
		wmat[wmat<0.001]<-min(wmat[wmat>0.001],na.rm=T)
		# get solution
		results<-find_closest_matrix_weighted(D=D, x=x_obs, weight_mat=wmat)
		Pfit<-convert_sol_to_matrix(results$X, nrow(D))
		colnames(Pfit)<-rownames(Pfit)<-colnames(D)
		# get the tournament matrix
		Hfit<-(Pfit+1)/2
		# get the log likelihood from binomial
		ll_mat<-D*log(Hfit)+t(D)*log(t(Hfit))
		# sum the upper triangle
		sum_ll <- sum(ll_mat[upper.tri(ll_mat)])
		# see if its more probably than the current best fit
		if(sum_ll>max_ll){
			max_ll <- sum_ll
			Pfit_best <-Pfit
		}

	}
	return(list(P=Pfit_best))
}


#### called within the quadratic programming results ####
# get outcome vector
get_vector <- function(P, n){
	Pv <- c()
	for(i in 1:(n-1)){
		Pv <- c(Pv, P[i, (i+1):n])
	}
	return (Pv)
}
# convert upper triangular solution to matrix
convert_sol_to_matrix <- function(sol, n){
	m <- length(sol)
	M <- matrix(0, n, n)
	j <- 1
	for(i in 1:(n - 1)){
		M[i, (i + 1): n] <- sol[j:(j + (n- i - 1))]
		j <- j + (n- i)
	}
	return (M - t(M))
}
# get the constraints for zero sum
get_constraints <- function(x, n, m){
	B <- matrix(0, n, m)
	for (i in 1:n){
		if (i < n){
			li <- 1 + (i - 1) * n  - i * (i - 1)/ 2
			B[i, li : (li + (n - i - 1))] <- x[(i + 1): n]
		}
		if( i > 1){
			if(i >= 3){
				s <- c(0, seq(n-2, n - (i - 1)))
			}
			if (i == 2){
				s <- c(0)
			}
			v <- cumsum(s)
			B[i, (i-1) + v] <- -x[1:(i-1)]
		}
	}
	A <- matrix(0, m, m)
	diag(A) <- 1
	return (rbind(B, A, -A))
}

######################################
#### Integrating dynamics


# replicator equation, not necessarily zero sum
replicator_eq <- function(time, x, params){
	with(as.list(params), {
		x[x<0]<-0
		x<-x/sum(x)
		dxdt<-x*(P%*%x-rep(t(x)%*%P%*%x,length(x)))
		return(list(dxdt))
	})
}


### integrating the intertidal dynamics
integrate_tidal_dynamics<- function(x0=NULL,pars,maxtime=1000,nsteps=100){
	x0<-x0/sum(x0)
	times <- seq(0, maxtime, length=nsteps)
	out <- as.matrix(ode(x0, times, replicator_eq, pars,method="ode45"))
	colnames(out)<-c("time",colnames(pars$P))
	return(out)
}

