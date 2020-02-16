DirectSampler <- function(N,n)
  {#used to directly sample N times from Ising of dim n (use when N small, otherwise use Gibbs)
  set.seed(1)
mu <- 0.1*rnorm(n) #creates mean matrix
set.seed(1)
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
set.seed(NULL)
S <- sample(I,N,replace=T,prob=Y) #picks a row of X according to pmf
J=c(1:N)
X <- 0*matrix(c(1:n*N),nrow = N, ncol = n)
for (j in J)
{
  X[j,] <- Q[S[j],1:n] # take jth sample as sth row of Y
}
return(X)
}