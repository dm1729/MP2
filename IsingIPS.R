IsingIPS <- function(x,M,n,eps) #x sample first moment, M sample second moment, n Graph size, eps accuracy
{
  #INITIALIZATION Want a distribution, edge set, and estimates for first and second moments
  #We intiialize at independent model with mean given by sample mean
mu <- x #n by 1 vector
C <- IsingDesign(n) #2^n by n matrix of all patterns
N=c(1:2^n)
I=c(1:n)
Y=2^(-n)*rep(1,2^n)
for (m in N) #loop to set up pmf for each sign pattern
{
    Y[m] <- Y[m]*prod((1-((-1)^(C[m,]>0))*mu)) #Mistake in paper as reqiuires binary vals in {0,1} not {-1,1}!
}
Z<-Y #REMOVE FROM CODE AND REPLACE RETURN WITH Y
#removed
Xi <- diag(n) #Initialize sample second moment
E <- diag(n) #Initialize E hat PUTING ONES ON DIAGONAL, SHOULDNT AFFECT ANYTHING
Delta <- 0*matrix(rep(1,n^2),nrow=n,ncol=n)
J <- 0*diag(n)
I <- c(1:(n-1))
#M[lower.tri(M)] <- 0 #We only update one half of Xi so we do this. M symmetric so no issue
Counter=0 #ADDED COUNTER TO EXIT LOOP FOR DIAGNOSTIC CHECKS
while ( (max(abs(mu-x))>=eps || any(Xi<M-10e-15) || max(abs((Xi-M)[E>0])>= eps) ) && Counter<10000 ) #Until all three of these do not hold.
{
  Counter=Counter+1
  for (i in I) #Initialize E+
  {
    II=c((i+1):n) #look over i<j
    for (j in II)
    {
      if (M[i,j]>x[i]*x[j]) #Has the effect of only looking through E plus (equivalent to the for i,j, in E+ step)
      {
        Delta[i,j] <- Delta(i,j,x,Y,M) #See Delta.R
        J[i,j] <- CanonicalJ(i,j,Y) #See CanonicalJ.R
        if (Delta[i,j]+J[i,j] > 0){
          #Update p by (14)
          for (m in N) #NOT CREATING A VALID PMF (sorted, wsn't signing the M_ij term)
            #SINCE EACH UPDATE USES WHOLE PMG IN DEMOINATOR, EARLIER UPDATES MESS UP LATER
            #INSTEAD NEED TO STORE NEW ONE SEPERATELY THEN UPDATE AT END (replaced with Z then assign Y<-Z)
          {
            Z[m] <- Y[m]*0.25*(1+C[m,i]*x[i]+C[m,j]*x[j]+C[m,i]*C[m,j]*M[i,j])/(sum(Y[C[,i]==C[m,i] & C[,j]==C[m,j]])) # denominator marginal
          }
          Y <- Z
          E[i,j] <- 1 #Puts an edge between i and j in the adjacency matrix
        } else {
          #Solve Delta_{ij}(lambda^*) equation
          LambdaStar <- LambdaStar(i,j,x,Y,M,J[i,j]) #J[i,j] not required if I've done this right!
          #Update p by (17)
          for (m in N)
          {
            Z[m] <- Y[m]*(0.25*(1+C[m,i]*x[i]+C[m,j]*x[j]+C[m,i]*C[m,j]*M[i,j]) + 0.25*C[m,i]*C[m,j]*LambdaStar)/(sum(Y[C[,i]==C[m,i] & C[,j]==C[m,j]])) # as before but add 0.25*x_ix_jlambda^*
          }
          E[i,j] <- 0
          Y <- Z
        } #end else
      } #end E plus restriction
    }
  } #end cycles through pairs
  
  #UPDATE mu and xi from p (eventually this means we stop repeating the above ['repeat' in paper] by exiting while loop)
  mu <- FirstMoment(Y)
  Xi <- SecondMoment(Y)
}#End while
return(Y) #return(Y,E,mu,Xi) Can't return multiple arguments, look over this.
}