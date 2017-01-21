# Derivative Function 
# 
# Takes the numerical derivative
# of a given function
#
# 01/18/2017
# Steven Novak
#
#

# Differentiate function
#
# Iterative function to calculate a functions derivative
# using a centered approximation and Richardsone extrapolation
#
# FUN is the function to differentiate
# X is the value to evaluate the derivative at
# err is the desired numeric precision

d2 <- function(FUN, X, err)
{
	maxIter <- 50
	step <- X 
	nextRow <- vector("numeric", maxIter)
	
	for(i in (1:maxIter))
	{
		nextRow[1] <- d2_cent(FUN, X, step)
		
		# ---- Richardson extrapolation routine
		k <- 2
		factor <- 2
		while( k <= i)
		{
			factor <- factor * 2
			nextRow[k] <- (factor * nextRow[k-1] - lastRow[k - 1]) / (factor - 1)
			k <- k + 1
		}

		# Check For change in precision < err
		if(i > 1)
		{
			error <- nextRow[i] - lastRow[i-1]
			if(abs(error) <= err)
			{
				print(i)
				return (nextRow[i])
			}
		}
		lastRow <- nextRow
		step <- step / 2
	}	
}

# Differentiate function
#
# Iterative function to calculate a functions derivative
# using a centered approximation and Richardsone extrapolation
#
# FUN is the function to differentiate
# X is the value to evaluate the derivative at
# err is the desired numeric precision

d <- function(FUN, X, err)
{
	maxIter <- 50
	step <- X 
	nextRow <- vector("numeric", maxIter)
	
	for(i in (1:maxIter))
	{
		nextRow[1] <- d_cent(FUN, X, step)
		
		# ---- Richardson extrapolation routine
		k <- 2
		factor <- 2
		while( k <= i)
		{
			factor <- factor * 2
			nextRow[k] <- (factor * nextRow[k-1] - lastRow[k - 1]) / (factor - 1)
			k <- k + 1
		}

		# Check For change in precision < err
		if(i > 1)
		{
			error <- nextRow[i] - lastRow[i-1]
			if(abs(error) <= err)
			{
				print(i)
				return (nextRow[i])
			}
		}
		lastRow <- nextRow
		step <- step / 2
	}
	
}


# Forward Looking Numerical Derivative
# FUN is a function where Y = FUN(X)
# X in the input(s) to FUN(X)
# Step is the dX to use in the derivative
#
d_for <- function(FUN, X, Step)
{
	yPlus <- FUN(X + Step)
	y <- FUN(X)
	return ((yPlus - y) / Step)
}

# Backward Looking Numerical Derivative
# FUN is a function where Y = FUN(X)
# X in the input(s) to FUN(X)
# Step is the dX to use in the derivative
#
d_back <- function(FUN, X, Step)
{
	yMin <- FUN(X - Step)
	y <- FUN(X)
	return ((y - yMin) / Step)
}

# Centered Numerical Derivative
# FUN is a function where Y = FUN(X)
# X in the input(s) to FUN(X)
# Step is the dX to use in the derivative
#
d_cent <- function(FUN, X, Step)
{
	yPlus <- FUN(X + Step)
	yMinus <- FUN(X - Step)
	return ((yPlus - yMinus) / (2 * Step))
}

# Forward Onesided Numerical Derivative
# FUN is a function where Y = FUN(X)
# X in the input(s) to FUN(X)
# Step is the dX to use in the derivative
#
d_fs <- function(FUN, X, Step)
{
	y <- FUN(X)
	yPlus <- FUN(X + Step)
	yPlus2 <- FUN(X + 2 * Step)
	return ((4 * yPlus - yPlus2 - 3 * y) / (2 * Step))
}

# Backward Onesided Numerical Derivative
# FUN is a function where Y = FUN(X)
# X in the input(s) to FUN(X)
# Step is the dX to use in the derivative
#
d_bs <- function(FUN, X, Step)
{
	y <- FUN(X)
	yMinus <- FUN(X - Step)
	yMinus2 <- FUN(X - 2 * Step)
	return ((3 * y - 4 * yMinus + yMinus2) / (2 * Step))
}

# Centered Second Numerical Derivative
# FUN is a function where Y = FUN(X)
# X in the input(s) to FUN(X)
# Step is the dX to use in the derivative
#
d2_cent <- function(FUN, X, Step)
{
	y <- FUN(X)
	yMinus <- FUN(X - Step)
	yPlus <- FUN(X + Step)
	return ((yPlus - 2 * y + yMinus) / (Step * Step))
}