function <- FlatteningRank(P) #Input a p.m.f. P should be a column vector of size 2^n
n=round(log2(P)) #ensures integer
#CREATES ARRAY such that P[1,2,1]=p_{-1,+1,-1} etc
dim(P) <- rep(2,n)
A <- aperm(P,c(4,3,2,1))
#FLATTENINGS
