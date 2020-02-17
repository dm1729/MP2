function <- FlatteningRank(P) #Input a p.m.f. P should be a column vector of size 2^n
n=round(log2(P)) #ensures integer
I=c(2:n-1) #Need to check all subsets of sizes 2 up to n-1
for (i in I)
{
  #need to cycle through the n (choose) i subsets of this size. Only take ones with 1 in (rest done on other side)
  
  #SHELVED FOR NOW< see FlatteningRank4 for small case