##Generating data from latent factor model

##library
library(MASS)

generateData = function( n, p, k, rep ){

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

# Generate data
temp2 = matrix( rnorm( N * p), N, p )
dat = tcrossprod(temp2, chol(Ot))

ktr = k
rktr = qr(Lambda)$rank
Lamtr = Lambda

#return data as output
return(list("data" = dat, "Var" = Ot, "repitition" = rep, "n" = n, "p" =p,
            "k.train" = ktr, "rank.train" = rktr, "Lambda.train" = Lamtr))

}
