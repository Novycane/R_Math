# Matrix Math Routines
#
# Routines for solving matrix problems
#
# 01/18/2017
#	- Added CROUT routine for LU factorization 
#
# Steven Novak
#

# --------------------------------------------- Helper Functions To Split Matrix

# Returns the lower matrix
# LU is a compact form L and U matrix where U is
# an upper unity matrix
L <- function(LU)
{
	n <- dim(LU)[1]
	l <- matrix(0,n,n)

	for(i in (1:n))
	{
		for(j in (1:i))
		{
			l[i,j] = LU[i,j]
		}
	}
	return (l)
}

# Returns the upper unity matrix
# LU is a compact form L and U matrix where U is
# an upper unity matrix
#
U <- function(LU)
{
	n <- dim(LU)[1]
	u <- matrix(0,n,n)

	for(i in (1:n))
	{
		for(j in (i:n))
		{
			if(i == j)
			{
				u[i,j] = 1
			}
			else
			{
				u[i,j] = LU[i,j]
			}
		}		
	}
	return (u)
}

# --------------------------------------------- Numerical Routines

# Cholesky Factorization
# Factors a symmetric positive definite matrix A
# 
#
Cholesky <- function(A)
{
	n <- dim(A)[1]

	for(i in 1:(n-1))
	{
		A[i,i] <- sqrt(A[i,i])

		l <- 0
		for(j in (i+1):n)
		{
			A[j,i] <- A[j,i] / A[i,i]
			A[i,j] <- 0
		}

		for(j in (i+1):n)
		{
			for(k in (i+1):n)
			{
				A[j,k] <- A[j,k] - A[k,i] * A[j,i]
			}
		}
	}
	A[n,n] <- sqrt(A[n,n])
	return(A)
}


# LU Factorization
# Factors a square non-singular matrix into
# a lower (L) and unit upper (U) matrix
#
# Output is a square matrix with the lower
# half being the L matrix, including the diagonal
# The upper half is the remaining matrix, with
# diagonal assumed to be 1's
#
LU <- function (X)
{
	rows <- dim(X)[1]
	columns <- dim(X)[2]

	if(rows != columns)
	{
		return(0)
	}

	lu <- matrix(0,rows,columns)

	for(j in (1:rows))
	{
		for(i in (j:columns))
		{
			tempSum = 0
			for(k in (1:j))
			{
				if(k == j)
				{
					u = 1;
				}
				else
				{
					u = lu[k,j]
				}
				tempSum = tempSum + lu[i,k] * u
			}
			lu[i,j] = X[i,j] - tempSum
		}

		for(i in (j:columns))
		{
			tempSum = 0
			for(k in (1:j))
			{
				if(i == k)
				{
					u = 1;
				}
				else
				{
					u = lu[k,i]
				}
				tempSum = tempSum + lu[j,k] * u
			}

			if(i != j)
			{
				lu[j,i] = (X[j,i] - tempSum) / lu[j,j]

			}
		}
	}
	return (lu)
}



# Back substitution routine
#
#

BackSub <- function(L, b)
{
	n <- dim(L)[1]
	x <- matrix(0,n,1)

	for(i in (1:n))
	{
		c = L[i,]
	
		x[i] <- (b[i] - (c %*% x)) / L[i,i]
	}
	return (x)
}


# Forward substitution routine
#
#
ForwardSub <- function(U, b)
{
	n <- dim(U)[1]
	x <- matrix(0,n,1)

	for(i in (n:1))
	{
		c = U[i,]

		x[i] <- (b[i] - (c %*% x)) / U[i,i]
	}
	return (x)
}

# Solve a linear system of equations Ax = b
# using LU decomposition
#


SolveLU <- function(A, b)
{
	lu <- LU(A)

	l <- L(lu)
	u <- U(lu)

	y <- BackSub(l, b)
	return (ForwardSub(u, y))
}

# Invert a matrix using LU Decomposition
#
#

InvertLU <- function(A)
{
	n <- dim(A)[1]
	lu <- LU(A)
	l <- L(lu)
	u <- U(lu)

	I <- matrix(0,n,n)
	for(i in (1:n))
	{
		I[i,i] = 1

		b <- I[i,]
		y <- BackSub(l, b)
		z <- ForwardSub(u, y)

		I[i,] <- t(z)
	}

	return(I)
}

# Find the determinate of a matrix using
# LU decomposition
#

DetLU <- function(A)
{
	lu <- LU(A)
	prod <- 1
	for(i in 1:dim(A)[1])
	{
		prod = prod * lu[i,i]
	}
	return (prod)
}

# Find the most dominate eigen value / vector
# using the power method
#

EigenPower <- function(A, X)
{
	n <- dim(A)[1]
	x <-X

	lambdaOld <- 10000
	lambda <- 0
	for(i in 1 : 100)
	{
		x <- A %*% x

		lambdaOld2 <- lambdaOld
		lambda <- 0
		for(j in 1:n)
		{
			if(abs(x[j]) > abs(lambda))
			{
				lambda <- x[j]
			}
		}

		x <- x / lambda
		err <- abs((lambda - lambdaOld))
		if(err <= 1e-5)
		{
		       break
		}
	}
	eval.parent(substitute(	X <- x))

	return(lambda)
}
