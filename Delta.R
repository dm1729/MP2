Delta <- function(i,j,x,Y,M) #Takes in the indices of vertices, means, and current iteration for distribution.
{
  n <- round(length(x)) #gives n
  C <- IsingDesign(n)
  q <- 0*diag(2)
  I=c(1:2)
  for (k in I)
  {
    for (l in I)
    {
      a <- 2*k-3
      b <- 2*l-3
      q[k,l] <- 0.25*(1+a*x[i]+b*x[j]+a*b*M[i,j])/(sum(Y[C[,i]==a & C[,j]==b])) #cf update for p in IsingIPS.R Replaced x_i, x_j with a,b
    }
  }
      
  Delta <- 0.25*(log(q[2,2])+log(q[1,1])-log(q[2,1])-log(q[1,2])) #Gives Delta as in (15)
  return(Delta)
}