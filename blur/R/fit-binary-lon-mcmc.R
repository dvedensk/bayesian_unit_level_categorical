#' This function fits the PL-MB model via Gibbs sampling
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
#' @param timepoints list of timepoints
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors

fit_bin_lon_mcmc <- function(Y, X, Psi, weights, timepoints, 
                             n_iter=2000, n_burn=500,
                             init_vals=NULL, hyperparams=NULL){

  n_keep <- n_iter - n_burn
  r <- ncol(Psi)
  n_tp <- max(timepoints) - min(timepoints) + 1
  p <- ncol(X)
  nn <- nrow(X) #number of responses (including repeat respondents)

  I_p <- Diagonal(p)
  I_r <- Diagonal(r)
  zero_p <- rep(0, p)
  zero_r <- rep(0, r)

  default_init_vals <- list(beta = rep(0, p),
                            eta = matrix(0, nrow=r, ncol=n_tp),
                            sigma2_eta = 1,
                            sigma2_eta = 1,
                            sigma2_eta_1 = 1,
                            phi=.5)

  init_vals <- validate_list_arg(default_init_vals, init_vals)

  default_hyperparams <- list(sigma2_beta = 3,
                              a=.1,
                              b=.1)
  
  hyperparams <- validate_list_arg(default_hyperparams, hyperparams)

  eta_out <- array(NA, dim = c(n_iter, r, n_tp))
  betas_out <- matrix(NA, n_iter, p)
  omegas_out <- array(NA, dim=c(n_iter, nn))
  phi_out <- sigma2_eta_out <- sigma2_eta_1_out <-rep(NA, n_keep)

  sigma2_eta <- init_vals$sigma2_eta
  sigma2_eta_1 <- init_vals$sigma2_eta_1
  phi <- init_vals$phi
    
  if(!all.equal(eta, default_init_vals$eta)) {
    for(tt in 1:n_tp){
      eta[,tt] <- MASS::mvrnorm(1, rep(0, r), diag(rep(.1, r))) ##
    }
  }

  kappa <- weights*(Y - .5)

  pb <- utils::txtProgressBar(min=0, max=n_iter, style=3)
  for(i in 1:n_iter){
    Psi_eta <- (Psi %*% eta)[cbind(seq_along(timepoints), timepoints)]
  
    omega <- BayesLogit::rpg.gamma(num=nn, h=weights,
                       z=as.numeric(X%*%beta + Psi_eta))
    Omega <- Diagonal(x=omega)

    sigma2_eta_1 <- extraDistr::rinvgamma(n=1, 
					  alpha=hyperparams$a + r/2,
                                          beta=hyperparams$b + t(eta[,1])%*%eta[,1]/2)
    
    sigma2_eta_shape <- hyperparams$a + r*(n_tp-1)/2
    eta_diff <- eta[,2:n_tp] - phi * eta[,1:(n_tp-1)]
    sigma2_eta_scale <- hyperparams$b + sum(eta_diff * eta_diff)/2

    sigma2_eta <- extraDistr::rinvgamma(n=1, 
					alpha=sigma2_eta_shape,
                                        beta=sigma2_eta_scale)
 
    denom <- sum(eta[,1:(n_tp-1)] * eta[,1:(n_tp-1)]) 
    mu_phi <- sum(eta[,2:n_tp] * eta[,1:(n_tp-1)])/denom
    sigma_phi <- sqrt(sigma2_eta/denom)
    phi <- phi_out[i] <- extraDistr::rtnorm(n=1, a=-1, b=1, mean=mu_phi, sd=sigma_phi)
    
    sqrt_XOmega <- sqrt(Omega)%*%X
    prec_beta <- t(sqrt_XOmega)%*%(sqrt_XOmega) + I_p/hyperparams$sigma2_beta
    mu_beta <- t(X)%*%Omega%*%(kappa/omega-Psi_eta)
    mu_beta <- solve(prec_beta, mu_beta)
    beta <- beta_out[i, ] <- rmvnorm_prec(n=1, mean=mu_beta, prec=prec_beta)

    ids_1 <- which(timepoints==1)
    Psi_1 <- Psi[ids_1,]
    omega_1 <- omega[ids_1]
    Omega_1 <- Diagonal(x=omega_1)
    X_1 <- X[ids_1,]
    kappa_1 <- kappa[ids_1]

    prec_eta_1 <- t(Psi_1)%*%Omega_1%*%Psi_1 + 
                         (1/sigma2_eta_1 + phi^2/sigma2_eta)*I_r 

    mu_eta_1 <- t(Psi_1)%*%Omega_1%*%(kappa_1/omega_1 - X_1%*%beta) +
        phi/sigma2_eta*eta[,2]
    mu_eta_1 <- solve(prec_eta_1, mu_eta_1)
    
    eta[, 1] <- eta_out[i,, 1] <- rmvnorm_prec(n=1, mean=mu_eta_1, prec=prec_eta_1)

    for(tt in 2:(n_tp-1)){
      ids_t <- which(timepoints==tt)
      Psi_t <- Psi[ids_t,]
      omega_t <- omega[ids_t]
      Omega_t <- Diagonal(x=omega_t)
      X_t <- X[ids_t,]
      kappa_t <- kappa[ids_t]

      prec_eta_t <- t(Psi_t)%*%Omega_t%*%Psi_t +
          (1+phi^2)/sigma2_eta * I_r
      mu_eta_t <- t(Psi_t)%*%Omega_t%*%(kappa_t/omega_t - X_t%*%beta) +
          phi/sigma2_eta*(eta[,(tt-1)] + eta[,(tt+1)])
      mu_eta_t <- solve(prec_eta_t, mu_eta_t)
      
       eta[, tt] <- eta_out[i,, tt] <- rmvnorm_prec(n=1, mean=mu_eta_t, prec=prec_eta_t)
    }

    ids_T <- which(timepoints==n_tp)
    Psi_T <- Psi[ids_T,]
    omega_T <- omega[ids_T]
    Omega_T <- Diagonal(x=omega_T)
    X_T <- X[ids_T,]
    kappa_T <- kappa[ids_T]

    prec_eta_T <- t(Psi_T)%*%Omega_T%*%Psi_T + I_r/sigma2_eta
    mu_eta_T <- t(Psi_T)%*%Omega_T%*%(kappa_T/omega_T - X_T%*%beta) +
        phi/sigma2_eta*eta[,(n_tp-1)]
    mu_eta_T <- solve(prec_eta_T, mu_eta_T)
  
    eta[, n_tp] <- eta_out[i,, tt] <- rmvnorm_prec(n=1, mean=mu_eta_T, prec=prec_eta_T)
      
    utils::setTxtProgressBar(pb, i)
  }
  posteriors <- list(posteriors = list(eta=eta_out[-c(1:n_burn),,],
                                       betas=beta_out[-c(1:n_burn), ],
                                       phi=phi_out[-c(1:n_burn)],
                                       sigma2_eta=sigma2_eta_out[-c(1:n_burn)], 
                                       sigma2_eta_1=sigma2_eta_1_out[-c(1:n_burn)]),
                     hyperparams=hyperparams)
  return(posteriors)
}


