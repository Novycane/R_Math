# Derivative Function 
# 
# Takes the numerical derivative
# of a given function
#
# 01/18/2017
# Steven Novak
#
#

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
d2 <- function(FUN, X, Step)
{
	y <- FUN(X)
	yMinus <- FUN(X - Step)
	yPlus <- FUN(X + Step)
	return ((yPlus - 2 * y + yMinus) / (Step * Step))
}