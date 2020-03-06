SignedIsingIPS <- function(X,eps) #Input dataset X and precision. Output MLE by fiting the various blobs.
{
n <- length(X[1,])
A <- IsingDesign(n-1) #Gives sign patterns which will be used for elemntwise multiplication to flip signs (coulduse bincom)
N <- length(X[,1]) #Returns sample size
Y <- matrix(0,nrow=2^(n-1),ncol=2^n)
S <-matrix(0,nrow=2^(n-1),ncol=n)
for (i in c(1:2^(n-1)) ){
  S[i,] <- c(1,A[i,]) #Creates 'sign flipper' vector. Never flip first one as only need to do half of possible flips
  Z <- t(t(X)*S[i,]) #Multiplies by sign flipper
  Y[i,] <- IsingIPS(colMeans(Z),SampleSecondMoment(Z),n,eps) #This pmf isnt actually what we want, need to relabel (or relabel Counts)
  #gives pmf for transformed data, we want pmf for original data (CANCELLED)
  Counts <- IsingCounts(Z) #Gives the counts for the transformed data
  L[i] <- Counts%*%log(Y[i,]) #Instead of using SignedMLE, we don't transform the pmf
  #But instead transform the counts
  #This should basically just permute the terms of the sum compared with previous method
  #Problem is then we have to store the sign flipper matrix to undo it...
  #Although perhaps doesn't matter just when ivnestigating flattening rank?
}
j<-which.max(L)
return(list(S[j,],Y[j,]))
}# ISSUE WITH RELABELLING Y BEFORE PROCESSING