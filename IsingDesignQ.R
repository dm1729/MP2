IsingDesignQ <- function(n) #Makes Ising design matrix including quadratic effects
{
  Z=IsingDesign(n) #Imports 2^n by n matrix
  X=0*matrix(c(1:2^(n-1)*(n*(n+1))),nrow=2^n,ncol=((n*(n+1))/2)) #Prepares quad effects matrix
  X[,1:n]=Z
  I=c(1:((n*(n-1))/2))
  M=c(1:(n-1))
  for (m in M)
  {
    X[,(n+1+(m-1)*n-(m*(m-1)/2)):(n+(m*n-(m*(m+1)/2)))] <- X[,m]*X[,(m+1):n] #Adds interactions as X1X2,X1X3,..,X2X3,..,X(n-1)Xn
  }
  return(X)
}