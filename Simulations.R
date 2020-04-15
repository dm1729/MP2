#INCLUDES SIMULATIONS RUN FOR WRITE UP
library(IsingIPS) #my package
library(rje)

n <- 5
d <- 100
eps <- 1e-8
FR <- rep(0,d)

for (i in c(1:d)){
  set.seed(i)
  x <- rdirichlet(1,rep(1,2^n))[1,]
  FR[i] <- FlatteningRank(SignedIsingIPScounts(x,eps)[[2]])
}

FR