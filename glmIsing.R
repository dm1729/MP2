glmIsing <- function(n,N) #N accuracy (c.10000) n nodes (c.5)
{
  #Fits a linear relationship with quadratic interactions between log-likelihood of joint and inputs
source("GibbsSampler.R")
source("IsingCounts.R")
source("IsingDesign.R")
source("IsingDesignQ.R")

Y <- IsingCounts(n,N)
X <- IsingDesignQ(n)

#fits ising model. Note interaction terms will be twice the lambda_{ij} since we have symmetric matrix
fittedIsing <- {glm(formula = Y ~ X, family = poisson())}
return(fittedIsing)
}