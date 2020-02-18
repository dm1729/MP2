CanonicalJ <- function(i,j,Y) #Takes in the indices of vertices as well as the total mass function. Computes J_{ij} by (12)
{
  N <- c(1:round(length(Y)))
  n <- round(log2(length(Y)))
  C <- IsingDesign(n)
  x <- C[N[C[,i]==1 & C[,j]==-1][1],] #Gives a choice for x
  y <- x
  y[i] <- -x[i]
  y[j] <- -x[j] #Makes y agree with x except on i,j
  a <- Y[which(apply(C, 1, function(z) all.equal(z[1:n], pmax(x,y))) == "TRUE")]
  b <- Y[which(apply(C, 1, function(z) all.equal(z[1:n], pmin(x,y))) == "TRUE")]
  c <- Y[which(apply(C, 1, function(z) all.equal(z[1:n], x)) == "TRUE")]
  d <- Y[which(apply(C, 1, function(z) all.equal(z[1:n], y)) == "TRUE")]
  J <- 0.25*(log(a)+log(b)-log(c)-log(d))
  return(J)
}