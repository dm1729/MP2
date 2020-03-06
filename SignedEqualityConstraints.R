SignedEqualityConstraints <- function(n,N,M) #Fits ISING IPS and fits M different dists from M Gibbs samples of size N
{
  I = c(1:M)
  v <- c(1:N)
  w <- c(1:n)
  R <- 0*I
  SV <- matrix(0,nrow=M,ncol=12)
  P <- matrix(0,nrow=M,ncol=n)
  S <- matrix(0,nrow=M,ncol=n*(n-1)/2)
  Y <- matrix(0,nrow=M,ncol=2^n)
  for (i in I){ #Produces M Gibbs Samples, each of size N and using a different base parametrisation
    set.seed(i)
    mu <- 0.1*rnorm(n) #creates mean matrix
    set.seed(i)
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
        p <- exp(mu[j]+2*x[-j]%*%lambda[j,-j])/(exp(-mu[j]-2*x[-j]%*%lambda[j,-j])+exp(mu[j]+2*x[-j]%*%lambda[j,-j]))
        x[j]<-2*rbinom(1,1,p)-1
      }
      X[k,]=x
    }
    P[i,] <- SignedIsingIPS(X,10e-3)[[1]] #matrix of flippers
    Y[i,] <- SignedIsingIPS(X,10e-3)[[2]] #matrix of pmfs
    R[i] <- FlatteningRank(Y[i,])
  }
  return(list(P,Y,R)) # List containing first moment, second moment, fitted pmf, flat rank, and flattening sing vals
}