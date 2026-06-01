#' This function fits the nominal (unordered categorical) model via Gibbs sampling
#' 
#' @import Matrix
#' @import extraDistr
#' @param X The input covariate matrix, one row per sample unit
#' @param Psi The matrix of spatial basis functions, one row per sample unit
#' @param Y A one-hot matrix of survey responses
#' @param n_iter The number of iterations
#' @param n_burn The length of burn in
#' @param weights The vector of survey weights
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors

fit_nominal_cs_mcmc <- function(Y, X, Psi, n_iter, n_burn,
                             weights=NULL, init_vals, hyperparams) {
  n <- nrow(Y)
  N <- rep(1, n)
  
  K <- ncol(Y)
  Nstar <- cbind(1, N - t(apply(Y, 1, cumsum)))[,1:(K-1)]
  weights <- rep(weights, K-1)
  posteriors <- list()
  for(kk in 1:(K-1)){
    nn <- Nstar[,kk]
    keep <- which(nn > 0)
    Y_k <- Y[keep, kk]
    X_k <- X[keep, ]
    Psi_k <- Psi[keep, ]
    weights_k <- weights[keep]
    
    posteriors[[kk]] <- fit_binary_cs_mcmc(X=X_k, Psi=Psi_k, Y=Y_k,
                                           hyperparams=hyperparams,
                                           init_vals=init_vals,
                                           n_iter=n_iter, n_burn=n_burn,
                                           weights=weights_k)$posteriors
  }
  return(list(posteriors=posteriors,
              hyperparams=hyperparams))
}
