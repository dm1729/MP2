GibbsSampler <- function(n,N)
{
#Gibbs sampler for Ising model of dimension n
#Specify number of iterations N, symmetric off-diagonal matrix lambda, mean vector mu
#n<-100
#N<-4
set.seed(NULL)
mu <- 0.1*rnorm(n) #creates mean matrix
set.seed(NULL)
lambda <- 0.1*matrix(rnorm(n^2),nrow=n) #creates interaction matrix. Only upper (or lower) diagonal is relevant.
lambda[lower.tri(lambda,diag=TRUE)] <- 0 #creates a lower diagonal matrix with zero diagonal
lambda <- lambda+t(lambda) #creates the related symmetric matrix
#Initialize sampler
x <- rep(1,n)
v <- c(1:N)
w <- c(1:n)
X <- 0*matrix(c(1:n*N),nrow=N,ncol=n) #nxN matrix of zeroes to include samples
set.seed(NULL) #Seemed to be keeping the set seed for the binomial generator as well without doing this.
  for (k in v)
  {
    for (j in w)
    {
    x[j]<-sign(rbinom(1,1,(exp(mu[j]+2*x[-j]%*%lambda[j,-j])/(exp(-mu[j]-2*x[-j]%*%lambda[j,-j])+exp(mu[j]+2*x[-j]%*%lambda[j,-j]))))-0.5)
    }
    X[k,]=x
  }
return(X)
}