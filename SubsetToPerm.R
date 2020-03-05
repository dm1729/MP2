SubsetToPerm <- function(s,n) #Input a subset of {1,...,n} Outputs corresponding perm.
{
  p <- rep(0,n)
i <- length(s)
p[1:i]<-s
for (j in c(1:n)){
  if (any(p==j)==FALSE){ #Not in there at all
    i <- i+1
    p[i] <- j #Add in j
  }
}
return(p)
}