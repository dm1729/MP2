#To be used as design matrix when counts are generates from IsingCounts.
IsingDesign <- function(n) #creates 2^n by n matrix
{
  X <- 0*matrix(c(1:2^n*n),nrow=2^n,ncol=n) #need a matrix 2^n by n
  I <- c(1:2^n)
  #v <- exp(c(1:n)*log(2))
  #M <- c(1:n)
  for (i in I){
      j <- i-1 #Initialize while loop
      while (j > 0) {
      a <- floor(log2(j))
        X[i,n-a] <- 1
      j <- j-2^a
      }
      

  }
  X[X==0] <- -1 #Replaes zeros with 1
  return(X)
} # E1071, function bin combinations