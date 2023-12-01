##Generating data from latent factor model
#' Generates data-sets from latent factor model
#'
#' @import MASS
#' @importFrom stats rnorm
#'
#' @param n Numeric, number of data points in each repetition
#' @param p Numeric, dimension of each data point
#' @param k Numeric, number of factors to generate from
#' @param rep Numeric, number of replicates in data-set
#'
#' @return list, containing data-sets generated from latent factor model with IID error with variance 0.01, true covariance matrix, Lambda and its rank.
#' The data replicates are ammended row-wise.
#'
#' @export
#'
#' @examples
#' newData = generateData( n = 1000, p = 50, k = 10, rep = 5 )

generateData = function( n, p, k, rep ){

  N = n * rep

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
  return(list("data" = dat, "Var" = Ot, "repilicate
              " = rep, "n" = n, "p" =p,
              "k.train" = ktr, "rank.train" = rktr, "Lambda.train" = Lamtr))
}
