# Numerical Integration Routines
#
# 01/19/2017
#
# Steven Novak
#

# --------------------------------------------- Numerical Routines

# Rhomberg Integration
#
# Implementation of Rhomberg integration routine
# using the CTrap routine below and a Richardson
# extrapolation routine to speed convergence
#

Rhomberg <- function(FUN, a, b, err)
{
	maxIter = 15
	h = b - a

	sum <- FUN(b) + FUN(a)
	lastRow <- vector("numeric", maxIter)
	nextRow <- vector("numeric", maxIter)
	points <- a

	for(j in (1:maxIter))
	{
		# ---- Trapezoidal integration routine
		h = h / 2
		len <- length(points)
		for(i in (1:len))
		{
			points <- c(points, points[i] + h)
			sum <- sum + 2 * FUN(points[i] + h)
		}
		nextRow[1] <- sum * h / 2

		# ---- Richardson extrapolation routine
		k <- 2
		factor <- 2
		while( k <= j)
		{
			factor <- factor * 2
			nextRow[k] <- (factor * nextRow[k-1] - lastRow[k - 1]) / (factor - 1)
			k <- k + 1
		}

		# Check For change in precision < err
		if(j > 1)
		{
			error <- nextRow[j] - lastRow[j-1]
			if(abs(error) <= err)
			{
				return (nextRow[j])
			}
		}
		lastRow <- nextRow
	}
}

# Composite Trapezoidal Integration
#
# Eatch iteration subdivides the interval a->b
# by 2. Only the new values from the subdivision
# are calculated for each iteration
#

CTrap <- function(FUN, a, b, err)
{
	maxIter = 15
	h = b - a

	sum <- FUN(b) + FUN(a)

	points <- a

	for(j in (1:maxIter))
	{
		h = h / 2
		len <- length(points)
		for(i in (1:len))
		{
			points <- c(points, points[i] + h)
			sum <- sum + 2 * FUN(points[i] + h)
		}
	}
}
