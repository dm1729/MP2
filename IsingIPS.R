IsingIPS <- function(x,M,n,eps) #x sample first moment, M sample second moment, n Graph size, eps accuracy
  #INITIALIZATION Want a distribution, edge set, and estimates for first and second moments
  #We intiialize at independent model with mean given by sample mean
mu <- x #n by 1 vector
C <- IsingDesign(n) #2^n by n matrix of all patterns
N=c(1:2^n)
I=c(1:n)
Y=2^(-n)*rep(1,2^n)
for (m in N) #loop to set up pmf for each sign pattern
{
  for (i in I)
  {
    Y[m] <- Y[m]*(1-((-1)^(C[m,i]>0))*mu[i]) #Mistake in paper as reqiuires binary vals in {0,1} not {-1,1}!
  }
}
Xi <- diag(n) #Initialize sample second moment
E <- 0*matrix(rep(1,n^2),nrow=n,ncol=n) #Initialize E hat
Delta <- 0*matrix(rep(1,n^2),nrow=n,ncol=n)
J <- 0*matrix(rep(1,n^2),nrow=n,ncol=n)
while max(abs(mu-x))>=eps | sum(Xi>=M)<2^n | max(abs(Xi-M))>= eps #Until all three of these do not hold.
{

for (i in I) #Initialize E+
{
  J=c((i+1):n) #look over i<j
  for (j in I)
  {
    if (M[i,j]>x[i]*x[j]) #Has the effect of only looking through E plus (equivalent to the for i,j, in E+ step)
    {
      Delta[i,j] <- Delta(i,j,x,Y,M) #See Delta.R
      J[i,j] <- CanonicalJ(i,j,Y) #See CanonicalJ.R
      if (Delta[i,j]+J[i,j] > 0)
      {
        #Update p by (14)
        for (m in N)
        {
          Y[m] <- Y[m]*0.25*(1+C[m,i]*x[i]+C[m,j]*x[j]+M[i,j])/(sum(Y[C[,i]==C[m,i] & C[,j]==C[m,j]])) # denominator marginal
        }
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
  mu <-
  Xi <-
}#End while

return(Y,E,mu,Xi)