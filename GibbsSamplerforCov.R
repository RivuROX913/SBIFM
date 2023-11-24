# Gibbs sampler for covariance estimation
# using mgps prior on factor loadings

library(Matrix)
# Define global constants
nrun <- 20000
burn <- 5000
thin <- 1
sp <- (nrun - burn) / thin  # Number of posterior samples

kinit <- p  # Number of factors to start with
b0 <- 1
b1 <- 0.0005
epsilon <- 1e-3  # Threshold limit
prop <- 1.00  # Proportion of redundant elements within columns

# Initialize arrays for storing results
mserep <- matrix(0, nrow = rep, ncol = 3)  # MSE, absolute bias (avg and max) in estimating cov matrix
mse1rep <- matrix(0, nrow = rep, ncol = 3)  # Same as above in the original scale in estimating cov matrix
nofrep <- matrix(0, nrow = rep, ncol = sp)  # Evolution of factors across replicates

for (g in 1:rep) {
  cat(paste("start replicate", g, "\n"))
  cat("--------------------\n")

  # Read data
  Y <- dat[((g - 1) * n + 1):(g * n), ]  # n x p data matrix
  M <- colMeans(Y)
  tY = t(Y) - M
  VY <- rowMeans(tY^2)
  Y <- t( tY / sqrt(VY) )
  #vhat = sum((t(Y) * sqrt(VY))^2)

  Ot1 <- Ot / as.numeric(sqrt(crossprod(VY)))  # True dispersion matrix of the transformed data
  num <- 0
  k <- kinit  # Number of factors to start with

  # Define hyperparameter values
  as <- 1
  bs <- 0.3
  df <- 3
  ad1 <- 2.1
  bd1 <- 1
  ad2 <- 3.1
  bd2 <- 1
  adf <- 1
  bdf <- 1

  # Initial values
  ps <- rgamma(p, as, bs)
  Sigma <- diag(1 / ps)  # Sigma = diagonal residual covariance
  Lambda <- matrix(0, nrow = p, ncol = k)
  eta <- matrix(rnorm(n * k), n, k)  # Factor loadings & latent factors
  meta <- matrix(0, nrow = n, ncol = k)
  veta <- diag(k)  # Latent factor distribution = standard normal
  psijh <- rgamma(p * k, shape = df / 2, scale = 2 / df)
  delta <- c(rgamma(1, shape = ad1, rate = bd1), rgamma(k - 1, shape = ad2, rate = bd2))
  tauh <- cumprod(delta)  # Global shrinkage coefficients
  Plam <- matrix(0, nrow = p, ncol = k)

  # Define output files specific to replicate
  nofout <- numeric(nrun + 1)  # Number of factors across iterations
  nofout[1] <- k
  nof1out <- numeric(sp)
  mseout <- matrix(0, nrow = sp, ncol = 3)  # Within a replicate, stores mse across MCMC iterations
  mse1out <- matrix(0, nrow = sp, ncol = 3)  # MSE in the original scale
  Omegaout <- rep(0, p^2)
  Omega1out <- rep(0, p^2)

  # Start Gibbs sampling
  for (i in 1:nrun) {
    # Update eta
    Lmsg <- Lambda * ps
    Veta1 <- diag(1, k) + t(Lambda) %*% Lmsg
    S = chol(Veta1)
    Veta <- chol2inv( S )
    Meta <- Y %*% Lmsg %*% Veta
    eta <- Meta + matrix(rnorm(n * k), n, k) %*% S  # Update eta in a block

    # Update Lambda (rue & held)
    eta2 <- t(eta) %*% eta
    for (j in 1:p) {
      Qlam <- diag(Plam[j, ]) + ps[j] * eta2
      blam <- ps[j] * t(eta) %*% Y[, j]
      Llam <- chol(Qlam)
      zlam <- matrix(rnorm(k), k, 1)
      vlam <- chol2inv(Llam)
      mlam <- vlam %*% blam
      ylam <- vlam %*% zlam
      Lambda[j, ] <- (ylam + mlam)
    }

    # Update psi_{jh}'s
    df2 = df / 2 + t(Lambda^2) * tauh / 2
    psijh <- matrix( 0, p , k  )
    for(j in 1:p)
    for(h in 1:k)
    psijh[j, h] = rgamma(1, df / 2 + 0.5, df2[h, j])

    # Update delta & tauh
    tmp <- colSums(Lambda^2 * psijh)
    ad <- ad1 + 0.5 * p * k
    bd <- bd1 + 0.5 * (1 / delta[1]) * crossprod(tauh ,tmp)
    delta[1] <- rgamma(1, ad, bd)
    tauh <- cumprod(delta)

    for (h in 2:k) {
      ad <- ad2 + 0.5 * p * (k - h + 1)
      bd <- bd2 + 0.5 * (1 / delta[h]) * crossprod(tauh[h:k] ,tmp[h:k])
      delta[h] <- rgamma(1, ad, bd)
      tauh <- cumprod(delta)
    }

    # Update Sigma
    Ytil <- Y - eta %*% t(Lambda)
    tmp2 = colSums(Ytil^2)
    for (j in 1:p)
    ps[j] <- rgamma(1, as + 0.5 * n, bs + 0.5 * tmp2[j])
    Sigma <- diag(1 / ps)

    # Update precision parameters
    Plam <- t( t(psijh) * tauh )

    # Make adaptations
    tmp3 = Ytil + Y
    S_h = rep(0, k)
    for(h in 1:k)
      S_h[h] = eta[,h] %*% tmp3 %*% Lambda[,h]
    mindex = which.min(S_h)
    num = sum(tmp2) + sum(S_h[-mindex])
    den = sum(Y^2)
    if( num / den > 0.999 )
    {
      k = k - 1
      Lambda = Lambda[, -mindex]
      eta = eta[, -mindex]
      delta = delta[-mindex]
      tauh <- cumprod(delta)  # Global shrinkage coefficients
      Plam <- Plam[, -mindex]
    }

    nofout[i + 1] <- k

    # Save sampled values (after thinning)
    if (i %% thin == 0 && i > burn) {
      Omega <- Lambda %*% t(Lambda) + Sigma
      Omega1 <- Omega * as.numeric(sqrt(crossprod(VY)))
      Omegaout <- Omegaout + as.vector(Omega) / sp
      Omega1out <- Omega1out + as.vector(Omega1) / sp
      nof1out[(i - burn) / thin] <- nofout[(i - burn) / thin] - num * (num > 0)
    }

    if (i %% 1000 == 0) {
      cat(i, "\n")
    }
  }

  # Summary measures specific to replicate
  # 1. Covariance matrix estimation
  errcov <- Omegaout - as.vector(Ot1)
  err1cov <- Omega1out - as.vector(Ot)
  mserep[g, ] <- c(mean(errcov^2), mean(abs(errcov)), max(abs(errcov)))
  mse1rep[g, ] <- c(mean(err1cov^2), mean(abs(err1cov)), max(abs(err1cov)))

  # 2. Evolution of factors
  nofrep[g, ] <- nof1out

  cat(paste("end replicate", g, "\n"))
  cat("--------------------\n")
}

# Display elapsed time

