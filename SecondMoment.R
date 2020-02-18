SecondMoment <- function(Y) #Input pmf, output second moment
{
  n <- round(log2(length(Y)))
  C <- IsingDesign(n)
  I <- c(1:n)
  Xi=0*diag(n)
  for (i in I)
    for (j in I)
      Xi[i,j] <- sum(Y[C[,i]*C[,j]==1])-sum(Y[C[,i]*C[,j]==-1]) #X_iX_j takes vals 1 or -1.
  return(Xi)
}