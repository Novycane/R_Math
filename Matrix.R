# Matrix Math Routines
#
# Routines for solving matrix problems
#
# 01/18/2017
#	- Added CROUT routine for LU factorization 
#
# Steven Novak
#

# --------------------------------------------- Numerical Routines

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