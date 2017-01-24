# Curve Fitting Routines
#
# Routines for curve fitting / regression
#
# 01/18/2017
#
# Steven Novak
#

# *								*
# *-------------------------------------------------------------*
# * Note-Function Should Be In The Form FUN(c(X, FittingParams))*
# *-------------------------------------------------------------*
# *								*

# --------------------------------------------- Dependencies
source("./Derivative.R")
source("./Matrix.R")

# --------------------------------------------- Numerical Routines


# Gauss-Newton
#
#
#

GN <- function(FUN, x, y, guess)
{
	maxIter <- 50
	
	for(iter in 1:maxIter)
	{
		oldGuess <- guess
		J <- Jacobian(FUN, x, guess, 1e-10)
		A <- t(J) %*% J
		A <- InvertLU(A)
		r <- Residual(f, x, y, guess)

		update <- c(0,A %*% (t(J) %*% r))
		guess <- guess + update
		#print(oldGuess - guess)
	}
	return(guess)
}

# Gradient Descent
#
#
#

GD <- function(FUN, x, y, guess)
{
	alpha <- 0.0001
	maxIter <- 20

	for(i in 1 : maxIter)
	{
		J <- Jacobian(FUN, x, guess, 1e-10)
		r <- Residual(FUN, x, y, guess)
		result <- c(0,t(J) %*% (r*r))
		guess <- guess + result
	}
	return(guess)
}


# Levenberg-Marquardt
#
#
#

LM <- function(FUN, x, y, guess)
{
	maxIter <- 50
	
	for(iter in 1:maxIter)
	{		
		H <- Hessian(FUN, guess, 1e-10)
		G <- Gradient(FUN, guess, 1e-10)
		A <- InvertLU(H)
		update <- A %*% G

		for(i in 1:dim(update)[1])
		{
			guess[i+1] <- guess[i+1] - update[i,1]
		}
		print(guess)
	}
	return(guess)	
}

# Jacobian for Gauss-Newton Routine
#
#
#

Jacobian <- function(FUN, x, params, err)
{
	n <- length(params) - 1
	m <- length(x)
	J <- matrix(0, m, n)
	
	for(i in 1:m)
	{
		params[1] <- x[i]
		for(j in 1:n)
		{
			val <- pd(FUN, params, j + 1, err)
			J[i, j] <- val  
		}
	}
	return(J)
}

# Hessian
#
#
#

Hessian <- function(FUN, params, err)
{
	n <- length(params) - 1
	J <- matrix(0,n,n)
	for(i in 1:n)
	{
		for(j in 1:n)
		{
			if(i == j)
			{
				J[i, j] <- pd2(FUN, params, i + 1, err)  
			}
			else
			{
				J[i, j] <- pd(FUN, params, i + 1, err) * pd(FUN, params, j + 1, err) 			
			}
		}
	}
	return(J)
}

# Gradient
#
#
#

Gradient <- function(FUN, params, err)
{
	n <- length(params) - 1
	G <- matrix(0,n,n)
	for(i in 1:n)
	{
		for(j in 1:n)
		{
			G[i, j] <- pd(FUN, params, i + 1, err)  
		}
	}
	return(G)
}

# Residual
#
#
#

Residual <- function(FUN, x, y, params)
{
	n <- length(x)
	r <- vector("numeric", n)
	for(i in 1:n)
	{
		params[1] <- x[i]
		r[i] <- y[i] - FUN(params)
	}
	return(r)
}