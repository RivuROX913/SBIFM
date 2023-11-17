##Generating data from latent factor model

library(MASS)

# Set parameters
rep <- 1
n <- 200
N <- n * rep
p <- 30
k <- 5

# Initialize variables
Lambda <- matrix(0, p, k)
numeff <- sample(k + 1:(2 * k), k, replace = FALSE)

# Generate loadings
for (h in 1:k) {
  temp <- sample(p, numeff[h], replace = FALSE)
  Lambda[temp, h] <- rnorm(numeff[h])
}

# Generative model: N(0, Lambda Lambda' + sigma^2 I)
mu <- rep(0, p)
Ot <- Lambda %*% t(Lambda) + 0.01 * diag(p)

# Generate data
dat <- mvrnorm(N, mu, Ot)

ktr <- k
rktr <- qr(Lambda)$rank
Lamtr <- Lambda

# Save data to a file
save(dat, Ot, rep, n, p, ktr, rktr, Lamtr, file = paste0("DataGen_p[", p, "]ktr[", k, "]rep[", rep, "]"))



