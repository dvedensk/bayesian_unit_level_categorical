#' Fit cross-sectional ordinal model with MCMC 
#'
#' @param Y n x K dimensional matrix of responses in one-hot format
#' @param X n x p design matrix of covariates
#' @param Psi n x r matrix of basis functions
#' @param weights n-dimensional vector of weights
#' @param n_iter number of iterations
#' @param n_burn number of iterations to discard for burn-in
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors
#' @returns a list of matrices of MCMC chains for each parameter

fit_ordinal_cs_mcmc <- function(Y,
                                X,
                                Psi,
                                weights = NULL,
                                n_iter=2000,
                                n_burn=500,
                                beta = NULL,
                                init_vals=NULL,
                                hyperparams=NULL) {

  if(!is.numeric(Y)) { stop("Y must be numeric") }
  if(!all(Y %in% c(0,1))) { stop("Y must be a vector of 0 and 1") }

  nn <- nrow(Y)
  N <- rep(1, nn)

  categories <- apply(Y, 1, which.max)

  p <- ncol(X)
  K <- ncol(Y)
  r <- ncol(Psi)

  I_p <- Diagonal(p)
  zero_p <- rep(0, p)
  I_r <- Diagonal(r)
  zero_r <- rep(0, r)
  I_k <- Diagonal(K-1)

  default_init_vals <- list(beta = rep(0, p),
                            eta = rep(0, r),
                            gamma = rep(0, K-1))

  init_vals <- validate_list_arg(default_init_vals, init_vals)

  default_hyperparams <- list(sigma2_beta = 3,
                              sigma2_eta = 1,
                              sigma2_gamma = 3,
                              a=.1,
                              b=.1)
  
  hyperparams <- validate_list_arg(default_hyperparams, hyperparams)

  if(is.null(weights)) weights <- rep(1,n)

  ###
  ### Setup storage
  ###
  n_keep <- n_iter - n_burn
  beta_post <- matrix(NA, n_keep, p)
  eta_post <- matrix(NA, n_keep, r)
  gamma_post <- matrix(NA, n_keep, K-1)
  sigma2_eta_post <- rep(NA, n_keep)

  keep <- 1:nn
  G_ind <- matrix(0, nrow=nn, ncol=(K-1))
  G_ind[,1] <- 1
  kappa <- as.matrix(weights * (Y[,-K] - 0.5)) #column for each cat
  kappa_block <- kappa[ , 1]
  for(k in 2:(K-1)){
      #indices for stacked matrices
      id_k <- which(categories >= k)
      keep <- c(keep, id_k) 

      #matrix for picking out the indices of gamma
      tmp <- matrix(0, nrow=length(id_k), ncol=(K-1))
      tmp[,k] <- 1
      G_ind <- rbind(G_ind, tmp)
 
      kappa_block <- c(kappa_block, kappa[id_k, k])    
  }
  kappa <- kappa_block
  X <- X[keep, ]
  Psi <- Psi[keep, ]
  weights <- weights[keep]

  eta <- init_vals$eta
  Psi_eta <- Psi%*%eta
  X_beta <- X%*%init_vals$beta
  gamma_ind <- G_ind %*% init_vals$gamma

  pb <- txtProgressBar(min=0, max=n_iter, style=3)
  for(iter in 1:n_iter){
     linpred <- X_beta + Psi_eta
     omega <- rpg(num=length(keep), h=weights, 
                  z=as.numeric(gamma_ind - linpred))
     Omega <- Diagonal(n=length(keep), x=omega)
     kappa_omega <- kappa/omega
     sigma2_eta <- extraDistr::rinvgamma(n=1,
                                         alpha=hyperparams$a + r/2,
                                         beta=hyperparams$b + t(eta)%*%eta/2)

     tPsiOmega <- t(Psi) %*% Omega
     precEta <- tPsiOmega %*% Psi + 1/sigma2_eta * I_r
     muEta <- tPsiOmega %*% (gamma_ind - kappa_omega - X_beta)
     muEta <- solve(precEta, muEta)
     eta <- rmvnorm_prec(n=1, mean = muEta, prec = precEta)
     Psi_eta <- Psi%*%eta

     ##Calculate params for block sampler
     tXOmega <-t(X) %*% Omega 
     muBeta <- tXOmega %*% (gamma_ind - kappa_omega - Psi_eta)
     precBeta <- tXOmega %*% X + 1/hyperparams$sigma2_beta * I_p
    
     tGOmega <- t(G_ind) %*% Omega
     muGamma <- tGOmega %*% (kappa_omega + X_beta + Psi_eta)
     precGamma <- tGOmega %*% G_ind + I_k*1/hyperparams$sigma2_gamma

     ##Sample as a block
     muBlock <- c(as.vector(muBeta), as.vector(muGamma))
     precBlock <- bdiag(precBeta, precGamma)
     muBlock <- solve(precBlock, muBlock)
     params <- rmvnorm_prec(n=1, mean=muBlock, prec=precBlock)

     beta <- params[c(1:p)]
     gamma <- params[-c(1:p)]

     X_beta <- X%*%beta
     gamma_ind <- G_ind %*% gamma
            
     if(iter > n_burn){
       beta_post[(iter-n_burn), ] <- beta
       eta_post[(iter-n_burn), ] <- eta
       gamma_post[(iter-n_burn), ] <- gamma
       sigma2_eta_post[(iter-n_burn)] <- sigma2_eta
     }
  setTxtProgressBar(pb, iter)
  }
  close(pb)
    
  posteriors <- list(posteriors = list(beta=beta_post,
                                       eta=eta_post,
                                       gamma=gamma_post,
                                       sigma2_eta=sigma2_eta_post),
                     hyperparams=hyperparams)
  return(posteriors)
}
