SampleSecondMoment <- function(D)
{
N <- length(D[,1])
n <- length(D[1,])
I <- c(1:N)
M <- 0*diag(n)
for (i in I){
  M <- M+D[i,]%*%t(D[i,])
}
M <- (1/N)*M
return(M)
}