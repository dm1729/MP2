EqualityConstraints <- function(n=4,N,M) #Fits ISING IPS and fits M different dists from M Gibbs samples of size N
{
I = c(1:M)
v <- c(1:N)
w <- c(1:n)
R <- 0*I
for (i in I){ #Produces M Gibbs Samples, each of size N and using a different base parametrisation
    set.seed(i)
    mu <- 0.01*rnorm(n) #creates mean matrix
    set.seed(i)
    lambda <- 0.1*matrix(rnorm(n^2),nrow=n) #creates interaction matrix. Only upper (or lower) diagonal is relevant.
    lambda[lower.tri(lambda,diag=TRUE)] <- 0 #creates a lower diagonal matrix with zero diagonal
    lambda <- lambda+t(lambda) #creates the related symmetric matrix
    #Initialize sampler
    x <- rep(1,n)
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
    A <- 0*c(1:n)
    for (l in c(1:n)){
        A[l] <- mean(X[,l]) #Sample mean of the Gibbs sample
    }
    B <- SampleSecondMoment(X) #Sample variance of the Gibbs
    Y <- IsingIPS(A,B,n,10e-3)
    R[i] <- FlatteningRank4(Y)
}
return(R)

}