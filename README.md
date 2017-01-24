# R_Math
Custom Numerical Math Routines in R

Library Contents

Matrix.R -- Routines to manipulate and do matrix calculations
-- LU uses the CROUT method to calculate the LU decomposition of a matrix
-- L returns the lower diagonal matrix after LU decomposition
-- U returns the uper diagonal matrix after LU decomposition
-- BackSub peforms the back substituion routine
-- ForwardSub performs the forward substitution routine
-- SolveLU solves a system of linear equations using LU decomposition
-- InvertLU inverts a matrix using LU decomposition
-- DetLU calculates the determinate using LU decomposition
-- EigenPower calculates the largest eigenvalue of a matrix
   -- May get moved to a serperate file in the future

Derivative.R -- Routines to approximate first and second order derivatives
-- d2 calculates the centered second derivate of a function to an arbitrary precision
-- d calculates the centered first derivate of a function to an arbitrary precision
-- d_for calculates the forward looking derivative of a function given a step size
-- d_back calculates the backward looking derivative of a function given a step size
-- d_cent calculates the centered derivative of a function given a step size
-- d_fs calculates the forward looking one sided derivative of a function given a step size
-- d_bs calculates the backward looking one sided derivative of a function given a step size
-- d2_cent calculates the centered second derivative of a function given a step size
-- pd2 calcuates the second partial derivative of a function to an arbitrary precision
-- pd calcuates the first partial derivative of a function to an arbitrary precision
-- pD2_cent calculates the second partial derivative od a function given a step size
-- pD_cent calculates the first partial derivative od a function given a step size

Integral.R -- Routines to calculate numeric integrals
-- Rhomberg calculates the integral of a function to an arbitrary precision using
   richardson extrapolation to speed convergence
-- CTrap calculates the integral of a function using the composit trapezoidal rule

CureFitting.R -- Routines for non linear curve fitting / regression
-- GN is the Gauss-Newton routine for non-linear curve fitting
-- GD is the gradient descent method for curve fitting
-- LM is the Levenberg-Marquardt method for non-linear curve fitting
-- Jacobian is a routine to calculate the Jacobian matrix of a function using finite difference
-- Hessian calculates the Hessian matrix of a function using finite difference
-- Gradient calculates the gradiant matrix of a function using finite difference
-- Residual calculates the residual of a function FUN(x)  given X and Y 