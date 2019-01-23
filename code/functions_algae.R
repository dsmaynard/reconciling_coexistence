

## lotka volterra dynamics for integrating and visualizing
LV_dyn<- function(time, x, params){
	with(as.list(params), {
		x[x<0]<-0
		dxdt<-x*(r+A%*%x)
		return(list(dxdt))
	})
}

# implement quadratic programming for Lotka-Volterra, supplied with target r, A, and allowable tolerance
fit_qp_LV <- function(A,r,x_obs,tol=0.5){
	nspp <- ncol(A)
	# convert A into a vector and create the target parameter vectors by appending r
	Avec <- as.numeric(t(A))
	parvec <- c(Avec,r)
	# number of pars to fit
	npar <- length(parvec)
	# get the bounds
	Alow <- Avec-abs(Avec*tol) 
	Aupp <- Avec+abs(Avec*tol) 
	rlow <- r-abs(r*tol)
	rupp <- r+abs(r*tol)
	# constrain the diagonals to be negative, and some small value away from zero
	Aupp[as.logical(as.numeric(diag(nspp)))] <- min(-1e-4,max(diag(A)))/4
	# inequality constraints
	h <- c(Alow,rlow,-Aupp,-rupp)
	G <- rbind(diag(npar),-diag(npar))
	# equality constraints, assuming at equilibirum
	E <- cbind(as.matrix(bdiag(replicate(nspp,matrix(x_obs,nrow=1),simplify = F))),diag(nspp))
	f <- rep(0,nrow(E))
	# fit the qp model, returning a silent warning if it doesn't converge
	fit <- tryCatch(lsei(A=diag(npar),B=parvec,E=E,F=f,G=G,H=h),  error=function(e) e, warning=function(w) w)
	return(fit)
}



# pass a vector or tolerances and cycle through, using qp to estimate the params, returning only those tolerances that converged
find_tolerance_qp_LV <- function(r, A, x_obs, alive=rep(TRUE,length(r)), tol_seq){
	good_tols <- NULL
	for(i in tol_seq){
		fit <- fit_qp_LV(A=A[alive,alive],r=r[alive],x_obs=x_obs[alive],tol=i)
		if(!is(fit,"warning") & !is(fit,"error")){
				good_tols <- rbind(good_tols,data.frame(tol=i,norm=fit$solutionNorm,error=fit$IsError))
		}
	}
	if(is.null(good_tols)){
		print("No Solutions!")
	}
	else{ 
		if(any(good_tols$error) & !all(good_tols$err)){
			print("Have solution, but some with error")
		}		
		if(all(good_tols$error)){
			print("Have solution, but all with error")
		}
	}
	return(good_tols)
}
