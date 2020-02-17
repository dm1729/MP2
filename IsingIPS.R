IsingIPS <- function(x,M,n,eps) #x sample first moment, M sample second moment, n Graph size, eps accuracy
  #INITIALIZATION Want a distribution, edge set, and estimates for first and second moments
  #We intiialize at independent model with mean given by sample mean
mu <- x #n by 1 vector
C <- IsingCounts(n) #2^n by n matrix
I=c(1:2^n)
M=c(1:n)
Y=2^(-n)*rep(1,n)
for (i in I) #loop to set up pmf for each sign pattern
{
  for (m in M)
  {
    Y[i]=Y[i]*(1-(-1)^(C[i,m])*mu(m))
  }
}
Xi <- diag(n) #Initialize sample second moment
E <- 0*matrix(rep(1,n^2),nrow=n,ncol=n) #Initialize E hat
Delta <- 0*matrix(rep(1,n^2),nrow=n,ncol=n)
J <- 0*matrix(rep(1,n^2),nrow=n,ncol=n)
while #Estimates are more than eps away or Xi<M
{

for i in M #Initialize E+
{
  J=c((i+1):n)
  for (j in J)
  {
    if (M[i,j]>x[i]*x[j]) #Has the effect of only looking through E plus (equivalent to the for i,j, in E+ step)
    {
      Delta[i,j]
      J[i,j]
      if (Delta[i,j]+J[i,j] > 0)
      {
        #Update p by (14)
        E[i,j] <- 1 #Puts an edge between i and j in the adjacency matrix
      } #end if
      else
      {
        #Solve Delta_{ij}(lambda^*) equation
        #Update p by (17)
        E[i,j] <- 0
      } #end else
    } #end E plus restriction
  }
} #end cycles through pairs

  #UPDATE mu and xi from p (eventually this means we stop repeating the above ['repeat' in paper] by exiting while loop)
  
}#End while