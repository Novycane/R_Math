# FiniteDiff.R
#
# Routines For Calculating solutions to differential equations
#
# 01/25/2017
#
#
# Steven Novak
#

# --------------------------------------------- Euler Routines

# Eulers Method for first order differential equations
#
# FUN a function in the form dY/dx = FUN(x)
# a is start of range
# b is the end of the range
# step is the step size
#

Euler <- function(FUN, a, b, step, ic)
{
	y <- vector("numeric", (b-a) / step)

	i = 1
	y[i] <- ic
	
	for(x in seq(a + step,b,step))
	{
		i = i + 1
		y[i] <- y[i-1] + step * FUN(y[i-1])
	}

	return(y)
}

# Eulers Method for first order differential equations
# using trapezoidal integration
#
# FUN a function in the form dY/dx = FUN(x)
# a is start of range
# b is the end of the range
# step is the step size
#

EulerTrap <- function(FUN, a, b, step, ic)
{
	y <- vector("numeric", (b-a) / step)

	i = 1
	y[i] <- ic
	
	for(x in seq(a + step,b,step))
	{
		i = i + 1
		yprime <- FUN(y[i-1])
		y[i] <- y[i-1] + (step / 2) * (yprime + FUN(y[i-1] + step * yprime))
	}

	return(y)

}