FlatteningRank <- function(P){ #Input a p.m.f. P should be a column vector of size 2^n
library(Matrix)
library(combinat)  
n <- round(log2(length(P))) #ensures integer
#CREATES ARRAY such that P[1,2,1]=p_{-1,+1,-1} etc
dim(P) <- rep(2,n) #P[2,1,2,1] = Pr(X_4=1,x_3=-1...)=p_{-1,1,-1,1} (inverted)
#FLATTENINGS
R <- rep(0,2^{n-1}-n-1)  #(2^n-2n-2)/2
i <- 0
for ( j in c( 2,(n-2) ) ){#Need to check subsets of this size
  for ( k in (choose((n-1),(j-1))) ){ #Choose the ones that aren't '1'
  p <- SubsetToPerm(combn(c(1:n),j)[,k],n) #First k columns gives the subsets we want
  A <- aperm(P,p) #Always fix one index.
  dim(A) <- c(2^j,2^(n-j)) #Flattens
  i <- i+1
  R[i] <- rankMatrix(A)
  }
}
return(max(R))
}