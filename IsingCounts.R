IsingCounts <- function(n,N)
{
#finds the count data C by representing each bin as binary string
I <- c(1:2^(n)) #check!
C <- 0*I
J <- c(1:N)
v <- round(exp(c((n-1):0)*log(2))) #Have to round here to get exact integers else equality condition breaks later
Z <- GibbsSampler(n,N) #Used to generate the counts (can replace with DirectSampler)
Z[Z < 0] <- 0 #Used to quickly convert between binary strings and representations. SO 10101 binary (21) is (1,-1,1,-1,1)
#Note the design matrix will be 2^n x (n+n(n-1)/2) to include the features for each combination
#Note that if the two way interactions were not included this would correspond to independence (can also see this as then we have n coeffecients to determine a binary vector 2^n possibilities)
for (i in I)
{
  for (j in J)
  {
  C[i] <- C[i] + (c(Z[j,]%*%v)==(i-1))
  }
}
return(C)
}