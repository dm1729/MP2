FirstMoment <- function(Y) #Input PMF, output mean
{
n <- round(log2(length(Y)))
C <- IsingDesign(n)
m <- Y%*%C #Creates mean vector by dotting each of the n columns of C with Y in turn
return(m)
}