LCMSampler <- function(n,N,l,p,q) #n dimension of LCM, N no. of samples p prob of H=1, p[1:n] probs of X_i=1|H=1, q for H=-1
{
S<-matrix(0,nrow=N,ncol=n)
for (i in c(1:N)){
  H <- rbinom(1,1,l)
  if (H==1){
    for (j in c(1:n)){
      #SET SEED
      S[i,j] <- rbinom(1,1,p[j])
    }
    }
  if (H==0){
    for (j in c(1:n)){
      #SET SEED
      S[i,j] <- rbinom(1,1,q[j])
    }
  }
  }

S[S==0] <- -1
return(S)
}