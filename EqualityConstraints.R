EqualityConstraints <- function(n=4,N,M) #Fits ISING IPS and fits M different dists from M Gibbs samples of size N
{
I = c(1:M)
v <- c(1:N)
w <- c(1:n)
R <- 0*I
SV <- matrix(0,nrow=M,ncol=12)
Q <- matrix(0,nrow=M,ncol=n)
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
    A <- colMeans(X)
    Q[i,] <- A
    B <- SampleSecondMoment(X) #Sample variance of the Gibbs
    for (m in c(1:(n-1))) #form taken from . Writes 
    {
      S[i,(1+(m-1)*n-(m*(m-1)/2)):((m*n-(m*(m+1)/2)))] <- B[m,(m+1):n] #Adds interactions as X1X2,X1X3,..,X2X3,..,X(n-1)Xn
    }
    Y[i,] <- IsingIPS(A,B,n,10e-3)
    R[i] <- FlatteningRank4(Y[i,])
    SV[i,] <- FlatteningRank4Mat(Y[i,]) #singular vals (in list of 12, 4 for each of the 3 flattenings)
}
return(list(Q,S,Y,R,SV)) # List containing first moment, second moment, fitted pmf, flat rank, and flattening sing vals
}