
# Set parameters
rep = 5
n = 200
N = n * rep
p = 30
k = 10

# Initialize variables
Lambda = matrix(0, p, k)
numeff = sample( (k + 1) : min( (3 * k), p), k, replace = FALSE)

# Generate loadings
for (h in 1:k) {
  temp = sample(p, numeff[h], replace = FALSE)
  Lambda[temp, h] = rnorm(numeff[h])
}

# Generative model: N(0, Lambda Lambda' + sigma^2 I)
mu = rep(0, p)
Ot = tcrossprod(Lambda) + 0.01 * diag(p)

library(testthat)
expect_equal(Ot, crossprod(chol(Ot)))

library(microbenchmark)
microbenchmark({
  temp2 = matrix( rnorm( N * p), N, p )
  dat = tcrossprod(temp2, chol(Ot))
},
dat = mvrnorm(N, mu, Ot))

