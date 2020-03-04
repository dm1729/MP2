IsingCountGen <- function(n,N) #Like Gibbs Sampler, but returns the counts rather than a matrix. Originally IsingCounts.R
{
#finds the count data C by representing each bin as binary string
I <- c(1:2^(n)) #check!
C <- 0*I
J <- c(1:N)
V <- round(exp(c((n-1):0)*log(2))) #Have to round here to get exact integers else equality condition breaks later
#BELOW THIS COPIED OVER FROM EQUALITYCONSTRAINTS GIBBS SAMPLER
v <- c(1:N)
w <- c(1:n)
set.seed(1)
mu <- 0.1*rnorm(n) #creates mean matrix
set.seed(1)
lambda <- 0.1*matrix(rnorm(n^2),nrow=n) #creates interaction matrix. Only upper (or lower) diagonal is relevant.
lambda[lower.tri(lambda,diag=TRUE)] <- 0 #creates a lower diagonal matrix with zero diagonal
lambda <- lambda+t(lambda) #creates the related symmetric matrix
#Initialize sampler
x <- rep(1,n)
X <- matrix(0,nrow=N,ncol=n) #nxN matrix of zeroes to include samples
set.seed(NULL) #Seemed to be keeping the set seed for the binomial generator as well without doing this.
for (k in v)
{
  for (j in w)
  {
    set.seed(2^(u)*3^(v)) # For reproducibility purposes.
    p <- exp(mu[j]+2*x[-j]%*%lambda[j,-j])/(exp(-mu[j]-2*x[-j]%*%lambda[j,-j])+exp(mu[j]+2*x[-j]%*%lambda[j,-j]))
    x[j]<-2*rbinom(1,1,p)-1
  }
  X[k,]=x
}
#END OF COPY
X[X < 0] <- 0 #Used to quickly convert between binary strings and representations. SO 10101 binary (21) is (1,-1,1,-1,1)
#Note the design matrix will be 2^n x (n+n(n-1)/2) to include the features for each combination
#Note that if the two way interactions were not included this would correspond to independence (can also see this as then we have n coeffecients to determine a binary vector 2^n possibilities)
for (i in I)
{
  for (j in J)
  {
  C[i] <- C[i] + (c(X[j,]%*%V)==(i-1))
  }
}
return(C)
}