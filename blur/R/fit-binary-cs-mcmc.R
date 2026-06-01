#' This function fits the binary cross-sectional model with MCMC
#' 
#' @import Matrix
#' @import BayesLogit
#' @param X The input covariate matrix, one row per sample unit
#' @param Psi The matrix of spatial basis functions, one row per sample unit
#' @param Y A vector of binary survey responses
#' @param n_iter The number of iterations
#' @param n_burn The length of burn in
#' @param weights The vector of survey weights
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors

fit_binary_cs_mcmc <- function(Y, X, Psi, weights, n_iter=1000, n_burn=500,
                               init_vals=NULL, hyperparams=NULL){
  #validate input data (Y should be a vector of 0s and 1s
  if(!is.numeric(Y)) { stop("Y must be numeric") }
  if(!all(Y %in% c(0,1))) { stop("Y must be a vector of 0 and 1") }

  p <- ncol(X)
  r <- ncol(Psi)
  n <- length(Y)

  default_init_vals <- list(beta = rep(0, p),
                            eta = rep(0, r),
                            sigma2_eta = 1)

  init_vals <- validate_list_arg(default_init_vals, init_vals)

  default_hyperparams <- list(sigma2_beta = 3,
                              sigma2_eta = 1,
                              a=.1,
                              b=.1)
  
  hyperparams <- validate_list_arg(default_hyperparams, hyperparams)

  if(is.null(weights)) weights <- rep(1,n)
  Binv <- (1/hyperparams$sigma2_beta) * Diagonal(p)
  kappa <- weights * (Y - 0.5)
  w <- rep(1, n)
    
  beta_out <- matrix(NA, nrow=n_iter, ncol=p)
  eta_out <- matrix(NA, nrow=n_iter, ncol=r)
  sigma2_eta_out <- rep(NA, n_iter)
  sigma2_eta <- init_vals$sigma2_eta
  Psi_eta <- Psi %*% init_vals$eta
  X_beta <- X %*% init_vals$beta
    
  pb <- utils::txtProgressBar(min=0, max=n_iter, style=3)
  for(i in 1:n_iter){
    ## Sample fixed effects
    precBeta <- t(X) %*% Diagonal(length(w), w) %*% X + Binv
    muBeta <- t(X) %*% Diagonal(length(w), w) %*% (kappa/w - Psi_eta)
    muBeta <- solve(precBeta, muBeta)
    beta <- beta_out[i,]  <- rmvnorm_prec(n=1,
                                          mean=muBeta,
                                          prec=precBeta)
    X_beta <- X %*% beta  

    ## Sample random effects
    Einv <- (1/sigma2_eta) * Diagonal(r)
    precEta <- t(Psi) %*% Diagonal(length(w), w) %*% Psi + Einv
    muEta <- t(Psi) %*% Diagonal(length(w), w) %*% (kappa/w - X_beta)
    muEta <- solve(precEta, muEta)
    eta <- eta_out[i,] <- rmvnorm_prec(n=1,
                                       mean=muEta,
                                       prec=precEta)
    #sum to zero
#    eta <- eta_out[i,]  <- eta - mean(eta)
    Psi_eta <- Psi %*% eta

    ## Sample RE variance
    sigma2_eta <- sigma2_eta_out[i] <- 1/stats::rgamma(n=1, shape=hyperparams$b + r/2,
                                                (hyperparams$a + 0.5*t(eta)%*%(eta)))
    
    ## Sample latent PG variables
    w <- BayesLogit::rpg(n, weights, as.numeric(X_beta + Psi_eta))
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
    
  return(list(posteriors = list(beta=beta_out[-c(1:n_burn),],
                                eta=eta_out[-c(1:n_burn),],
                                sigma2_eta=sigma2_eta_out[-c(1:n_burn)]),
              hyperparams=hyperparams))        
}
