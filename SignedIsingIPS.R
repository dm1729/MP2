SignedIsingIPS <- function(X,eps) #Input dataset X and precision. Output MLE by fiting the various blobs.
{
n <- length(X[1,])
A <- IsingDesign(n-1) #Gives sign patterns which will be used for elemntwise multiplication to flip signs (coulduse bincom)
N <- length(X[,1]) #Returns sample size
Y <- matrix(0,nrow=2^(n-1),ncol=2^n)
for (i in c(1:2^(n-1)) ){
  S <- c(1,A[i,]) #Creates 'sign flipper' vector. Never flip first one as only need to do half of possible flips
  Z <- t(t(X)*S) #Multiplies by sign flipper
  J <- IsingIPS(colMeans(Z),SampleSecondMoment(Z),n,eps) #This pmf isnt actually what we want, need to relabel (or relabel Counts)
  #gives pmf for transformed data, we want pmf for original data
  dim(J)=rep(2,n)
  J <- aperm(J,c(n:1)) #Tensorising J
  #FIDDLE WITH J
  for (j in c(1:n))
    K
  
  J <- aperm(J,c(n:1)) #Inverting transformation
  
  dim(J)=2^n
  Y[i,] <- J
  }
Counts <- IsingCounts(X) #Required for computation of likelihood
P <- SignedMLE(Y,Counts) #Gives the row of Y with the best fit in terms of likelihood.
return(P)
}# ISSUE WITH RELABELLING Y BEFORE PROCESSING