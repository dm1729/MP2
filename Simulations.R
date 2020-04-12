#INCLUDES SIMULATIONS RUN FOR WRITE UP
library(IsingIPS) #my package
library(rje) # for rdirichlet

GibbsSampler <- function(n,N)
{
  #Gibbs sampler for Ising model of dimension n
  #Specify number of iterations N, symmetric off-diagonal matrix lambda, mean vector mu
  #n<-100
  #N<-4
  set.seed(1)
  mu <- 0.1*rnorm(n) #creates mean matrix
  set.seed(1)
  lambda <- 0.1*matrix(rnorm(n^2),nrow=n) #creates interaction matrix. Only upper (or lower) diagonal is relevant.
  lambda[lower.tri(lambda,diag=TRUE)] <- 0 #creates a lower diagonal matrix with zero diagonal
  lambda <- lambda+t(lambda) #creates the related symmetric matrix
  #Initialize sampler
  x <- rep(1,n)
  v <- c(1:N)
  w <- c(1:n)
  X <- 0*matrix(c(1:n*N),nrow=N,ncol=n) #nxN matrix of zeroes to include samples
  #Seemed to be keeping the set seed for the binomial generator as well without doing this.
  for (k in v)
  {
    for (j in w)
    {
      set.seed(2^k*3^j)
      x[j]<-sign(rbinom(1,1,(exp(mu[j]+2*x[-j]%*%lambda[j,-j])/(exp(-mu[j]-2*x[-j]%*%lambda[j,-j])+exp(mu[j]+2*x[-j]%*%lambda[j,-j]))))-0.5)
    }
    X[k,]=x
  }
  return(X)
}

DirectSampler <- function(n,N)
{#used to directly sample N times from Ising of dim n (use when N small, otherwise use Gibbs)
  set.seed(2)
  mu <- 0.1*rnorm(n) #creates mean matrix
  set.seed(2)
  lambda <- 0.1*matrix(rnorm(n^2),nrow=n) #creates interaction matrix. Only upper (or lower) diagonal is relevant.
  lambda[lower.tri(lambda,diag=TRUE)] <- 0 #creates a lower diagonal matrix with zero diagonal
  lambda <- lambda+t(lambda) #creates the related symmetric matrix
  #Ising Design will give patterns. Then use likelihood to compute respective mass function
  Y <- 0*c(1:n)
  Q <- IsingDesignQ(n)
  I <- c(1:2^n)
  U <- 2*lambda[lower.tri(lambda)] #contains interactions in order X12, X13, X14, X23 etc... Times 2 to get symmetry
  for (i in I){
    Y[i] <- exp(Q[i,1:n]%*%mu + Q[i,(n+1):(n*(n+1)/2)]%*%U) #then need to renormalize (Maybe using IsingDesignQ and flatten lambda for interaction)
  }
  #Y[,n+1] <- (Y[,n+1]/sum(Y[,n+1])) #renormalizes to make a valid pmf (don't need to do actually)
  set.seed(2)
  S <- sample(I,N,replace=T,prob=Y) #picks a row of X according to pmf
  J=c(1:N)
  X <- 0*matrix(c(1:n*N),nrow = N, ncol = n)
  for (j in J)
  {
    X[j,] <- Q[S[j],1:n] # take jth sample as sth row of Y
  }
  return(X)
}

LCMSampler <- function(n,N,l,p,q) #n dimension of LCM, N no. of samples l prob of H=1, p[1:n] probs of X_i=1|H=1, q for H=-1
{
  S<-matrix(0,nrow=N,ncol=n)
  for (i in c(1:N)){
    set.seed(123+i)
    H <- rbinom(1,1,l)
    if (H==1){
      for (j in c(1:n)){
        #SET SEED
        set.seed(456+i+j)
        S[i,j] <- rbinom(1,1,p[j])
      }
    }
    if (H==0){
      for (j in c(1:n)){
        #SET SEED
        set.seed(789+i+j)
        S[i,j] <- rbinom(1,1,q[j])
      }
    }
  }
  
  S[S==0] <- -1
  return(S)
}

nX <- 4
nY <- 4
nZ <- 4

NX <- 100
NY <- 100
NZ <- 100

l <- 0.25
p <- rep(0.25,nZ)
q <- rep(0.35,nZ)

X <- GibbsSampler(nX,NX)
Y <- DirectSampler(nY,NY)
Z <- LCMSampler(nZ,NZ,l,p,q)

epsX <- 1e-4
epsY <- 1e-4
epsZ <- 1e-4

PX <- SignedIsingIPS(X,epsX)[[2]]
PY <- SignedIsingIPS(Y,epsY)[[2]]
PZ <- SignedIsingIPS(Z,epsZ)[[2]]
  
FX <- FlatteningRank(PX)
FY <- FlatteningRank(PY)
FZ <- FlatteningRank(PZ)
  
