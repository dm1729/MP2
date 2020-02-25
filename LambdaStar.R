LambdaStar <- function(i,j,x,Y,M,J=CanonicalJ(i,j,Y)) #Creates quadratic equation from delta_{ij}(lambda)=-J_{ij}
{
n <- round(length(x)) #gives n
C <- IsingDesign(n)
p <- 0*diag(2)
e <- p #THIS WAS A q!!!
I=c(1:2)
for (k in I)
{
  for (l in I)
  {
    a <- 2*k-3
    b <- 2*l-3
    e[k,l] <- 0.25*(1+a*x[i]+b*x[j]+a*b*M[i,j]) #e[1,1] is e_{ij}(-1,1) etc.
    p[k,l] <- (sum(Y[C[,i]==a & C[,j]==b])) #Sums over elements of the pmf with X_i=a and X_j=b
  }
}
R <- (p[1,1]*p[2,2])/(p[1,2]*p[2,1])*exp(-4*J)
a <- 1-R
b <- e[1,1]+e[2,2]+R*(e[1,2]+e[2,1])
c <- e[1,1]*e[2,2]-R*(e[1,2]*e[2,1])
if (a>0){
x <- (-b+sqrt(b^2-4*a*c))/(2*a)
lambda <- 4*x
}else{
  x <- (-c)/b
  lambda <- 4*x
}
return(lambda)
}