## Gibbs sampler for covariance estimation
# using mgps prior on factor loadings
#' Estimates loading matrix and covariance for factor model data through Sparse Bayesian Infinite Factor Model
#'
#' @import Matrix
#' @import Rcpp
#' @import RcppArmadillo
#' @include Datagen.R
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @useDynLib SBIFM, .registration = TRUE
#' @param data List, factor model data generated through generateData function
#' @param nrun Numeric, number of iterations to run the algorithm
#' @param burn Numeric, number of iterations to ignore
#' @param thin Numeric, thinning parameter for the outcome
#' @param epsilon proper fraction, used in truncation of factors
#'
#' @return List containing an estiate of the loading matrix, covariance matrix, error variance, error in estimating the covariance (if true covariance is supplied), number of chosen factors and posterior estimate of factors accross MCMC iterations.
#'
#' @references
#' Bhattacharya et al. (2011). "Sparse Bayesian Infinite Factor Models",
#' \emph{Biometrika},
#' <https://doi.org/10.1093/biomet/asr013>
#'
#' @references Schiavon, Lorenzo, and Antonio Canale. (2020).
#'  "On the truncation criteria in infinite factor models",
#'  \emph{Stat},
#'  <https://doi.org/10.1002/sta4.298>
#'
#' @details
#' The MCMC algorithm used is taken from Bhattacharya et. al (2011) with a truncation criteria derived from Schiavon et. al. (2020).
#'
#'
#' @export
#'
#' @examples
#' newdata = generateData( n = 100, p = 50, k = 10, rep = 2 )
#' GibbsCov(data = newdata, nrun = 2000, burn = 500, thin = 3, epsilon = 1e-4)


GibbsCov = function(data, nrun = 2000, burn = 500, thin = 3, epsilon = 1e-4)
{

  ##consistency checks
  stopifnot("nrun must be positive" = nrun > 0)
  stopifnot("burn must be positive" = burn > 0)
  stopifnot("nrun must be greater than burn" = nrun > burn)
  stopifnot("thin must be positive" = thin > 0)
  stopifnot("thin must be smaller than nrun - burn" = thin < nrun - burn)
  stopifnot("epsilon must be a proper proportion" = epsilon > 0 & epsilon < 1)

  ##definitions
  dat = data$data
  Ot = data$Var
  rep = data$replicate
  n = data$n
  p = data$p
  sp = (nrun - burn) / thin  # Number of posterior samples


  # Initialize arrays for storing results
  if(is.null(Ot) == FALSE)
  {
    mserep = matrix(0, nrow = rep, ncol = 3)  # MSE, absolute bias (avg and max) in estimating cov matrix
    mse1rep = matrix(0, nrow = rep, ncol = 3)  # Same as above in the original scale in estimating cov matrix
    colnames(mse1rep) = c("mean.sq", "mean.abs", "max.abs")
  }
  nofrep = matrix(0, nrow = rep, ncol = sp)  # Evolution of factors across replicates
  postfrep = rep(0, sp)
  Sigmarep = rep(0, p)

for (g in 1:rep) {
  cat(paste("start replicate", g, "\n"))
  cat("--------------------\n")

  # Read data
  Y = dat[((g - 1) * n + 1):(g * n), ]  # n x p data matrix
  M = colMeans(Y)
  tY = t(Y) - M
  VY = rowMeans(tY^2)
  Y = t( tY / sqrt(VY) )

  Ot1 = Ot / as.numeric(sqrt(crossprod(VY)))  # True dispersion matrix of the transformed data
  k = p  # Number of factors to start with

  # Define hyperparameter values
  as = 1
  bs = 0.3
  df = 3
  ad1 = 2.1
  bd1 = 1
  ad2 = 3.1
  bd2 = 1
  adf = 1
  bdf = 1

  # Initial values
  ps = rgamma(p, as, bs)
  Sigma = diag(1 / ps)  # Sigma = diagonal residual covariance
  Lambda = matrix(0, nrow = p, ncol = k)
  eta = matrix(rnorm(n * k), n, k)  # Factor loadings & latent factors
  meta = matrix(0, nrow = n, ncol = k)
  veta = diag(k)  # Latent factor distribution = standard normal
  psijh = rgamma(p * k, shape = df / 2, scale = 2 / df)
  delta = c(rgamma(1, shape = ad1, rate = bd1), rgamma(k - 1, shape = ad2, rate = bd2))
  tauh = cumprod(delta)  # Global shrinkage coefficients
  Plam = matrix(0, nrow = p, ncol = k)

  # Define output files specific to replicate
  nofout = numeric(nrun + 1)  # Number of factors across iterations
  nofout[1] = k
  nof1out = numeric(sp)
  if(is.null(Ot) == FALSE)
  {
    mseout = matrix(0, nrow = sp, ncol = 3)  # Within a replicate, stores mse across MCMC iterations
    mse1out = matrix(0, nrow = sp, ncol = 3)  # MSE in the original scale
  }
  Omegaout = Omega1out = matrix(0, p, p)

  # Start Gibbs sampling
  for (i in 1:nrun) {
    # Update eta
    eta = updateEta_c(Lambda, ps, k, Y, n)

    # Update Lambda
    Lambda = updateLambda_c(eta, Plam, ps, Y, k, p)

    # Update psi_{jh}'s
    psijh = updatePsi_c(df, Lambda, tauh, p, k)

    # Update delta & tauh
    out = updateDeltaTauh_c(Lambda, psijh, ad1, p, k, bd1, delta, tauh, ad2, bd2)
    delta = out$delta
    tauh = as.vector(out$tauh)

    # Update Sigma
    Ytil = Y - eta %*% t(Lambda)
    tmp2 = colSums(Ytil^2)
    ps = updateSigma_c(tmp2, p, as, n, bs)
    Sigma = diag(1 / as.vector(ps))

    # Update precision parameters
    Plam = t( t(psijh) * tauh )

    # Make adaptations
    if(k > 1) {
      tmp3 = Ytil + Y
      S_h = c()
      for(h in 1:k)
        S_h[h] = eta[,h] %*% tmp3 %*% Lambda[,h]
      mindex = which.min(S_h)
      num = sum(tmp2) + sum(S_h[-mindex])
      den = sum(Y^2)
      if(num / den > 1 - epsilon)
      {
        k = k - 1
        Lambda = Lambda[, -mindex]
        eta = eta[, -mindex]
        delta = delta[-mindex]
        tauh = cumprod(delta)  # Global shrinkage coefficients
        Plam = Plam[, -mindex]
      }
    }

    nofout[i + 1] = k

    # Save sampled values (after thinning)
    if (i %% thin == 0 && i > burn) {
      Omega = Lambda %*% t(Lambda) + Sigma
      Omega1 = Omega * as.numeric(sqrt(crossprod(VY)))
      Omegaout = Omegaout + Omega / sp
      Omega1out = Omega1out + Omega1 / sp
      nof1out[(i - burn) / thin] = nofout[(i - burn) / thin]
      Sigmarep = Sigmarep + as.numeric(sqrt(crossprod(VY))) * Sigma
    }
    if (i %% 1000 == 0) {
      cat(" " , i, "\n")
    }
  }

  if(is.null(Ot) == FALSE)
  {
    # 1. Error in covariance matrix estimation
    errcov = Omegaout - Ot1
    err1cov = Omega1out - Ot
    mserep[g, ] = c(mean(errcov^2), mean(abs(errcov)), max(abs(errcov)))
    mse1rep[g, ] = c(mean(err1cov^2), mean(abs(err1cov)), max(abs(err1cov)))
  }

  # 3. Covariance of errors
  Sigma1rep = Sigmarep / sp

  # 2. Evolution of factors
  nofrep[g, ] = nof1out

  # 3. posterior estimate of factors
  postfrep = colMeans(nofrep)

  cat(paste("end replicate", g, "\n"))
  cat("--------------------\n")
}
  if(is.null(Ot))
  {
    return(list( "Lambda" = Lambda, "Eta" = eta, "Cov" = Omega1out, "Sigma" = Sigma1rep,
                 "Factor" = nofrep[,sp-1], "post.factor" = postfrep[-sp]))
  } else {
    return(list( "Lambda" = Lambda, "Eta" = eta, "Cov" = Omega1out, "Sigma" = Sigma1rep,
                 "Error" = mse1rep, "Factor" = nofrep[,sp-1], "post.factor" = postfrep[-sp]))
  }
}
