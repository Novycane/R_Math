# Numerical Integration Routines
#
# 01/19/2017
#
# Steven Novak
#

# --------------------------------------------- Numerical Routines

# Rhomberg Integration
#
#
Rhomberg <- function(FUN, a, b, err)
{

}

# Composite Trapezoidal Integration
#
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
		print(sum * h / 2)
	}
	#print(points)
}