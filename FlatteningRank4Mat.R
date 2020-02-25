FlatteningRank4Mat <- function(P) #Input a PMF for a 4 dimensional binary model (16 masses)
{
library("Matrix", lib.loc="C:/Program Files/R/R-3.6.2/library")
I <- c(1:4)
# checking {1,2}
A <- 0*matrix(c(1:16),nrow=4,ncol=4)
for (i in I)
{
  for (j in I)
  {
    A[i,j] <- P[ 4*i-3 + j-1 ]
  }
}
# checking {1,3}
B <- 0*matrix(c(1:16),nrow=4,ncol=4)
for (i in I)
{
  for (j in I)
  {
    B[i,j] <- P[ 1+8*(i>(4/2))+2*((i/2)%%1==0) + (j-1)+2*(j>=3) ] #adds 8 when halfway through, adds 2 on even entries
  }
}
# checking {1,4} (actually think this is {2,3} if using convention of above but works out)
C <- 0*matrix(c(1:16),nrow=4,ncol=4)
for (i in I)
{
  for (j in I)
  {
    C[i,j] <- P[ 1+8*(i>(4/2))+((i/2)%%1==0) + 2*(j-1) ] # gives 1 2 9 10 etc
  }
}
r <- max(rankMatrix(A,tol=1e-8),rankMatrix(B,tol=1e-4),rankMatrix(C,tol=1e-4)) #finds maximum rank (built in tolerance seems like)
s=0*matrix(rep(1,12))#
s[1:4] <- svd(A)$d
s[5:8] <- svd(B)$d
s[9:12] <- svd(C)$d
return(s)
}
