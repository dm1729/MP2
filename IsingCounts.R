IsingCounts <- function(X) #Like isingCountGen but doesn't generate the raw output of e.g. GibbsSampler.R instead takes it in
{
  #ALL COPIED FROM BOTTOM OF ISINGCOUNTGEN
  n <- length(X[1,])
  C <- rep(0,2^n)
  V <- round(exp(c((n-1):0)*log(2)))
  X[X < 0] <- 0 #Used to quickly convert between binary strings and representations. SO 10101 binary (21) is (1,-1,1,-1,1)
  #Note the design matrix will be 2^n x (n+n(n-1)/2) to include the features for each combination
  #Note that if the two way interactions were not included this would correspond to independence (can also see this as then we have n coeffecients to determine a binary vector 2^n possibilities)
  for (i in c( 1:2^n ) )
  {
    for (j in c(1: length(X[,1])))
    {
      C[i] <- C[i] + (c(X[j,]%*%V)==(i-1))
    }
  }
  return(C)
}