SignedMLE <- function(Y,Counts) #Y 2^(n-1) * 2^(n) matrix of fitted pmf for each blob (one blob each row). Counts vector length 2^n
{
  n <- round(log2(length(Y[1,]))) # Returns n
L <- rep(0,2^(n-1)) # Allocates likelihood vector
  for (i in c(1:2^(n-1))){
    L[i] <- Counts%*%log(Y[i,])
  }
j <- which.max(L) # Returns index of L corresponding to max likelihood across blobs
return(Y[j,]) # Gives fitted pmf based on max multinomial likelihood.
}