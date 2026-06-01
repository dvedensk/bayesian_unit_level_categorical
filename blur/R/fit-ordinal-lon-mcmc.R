#' This function fits the longitudinal ordinal model via MCMC
#' 
#' @import Matrix
#' @import BayesLogit
#' @param X The input covariate matrix, one row per sample unit
#' @param Psi The matrix of spatial basis functions, one row per sample unit
#' @param timepoints A vector of the time at which measurements were taken
#' @param Y A vector of binary survey responses
#' @param n_iter The number of iterations
#' @param n_burn The length of burn in
#' @param weights The vector of survey weights
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors
#' @param prev_covar List of previous responses if longitudinal effect is desired

fit_ordinal_lon_mcmc <- function(Y,
                                 X,
                                 Psi,
#                                 areas,
                                 weights,
                                 timepoints,
                                 n_iter=2000,
                                 n_burn=500,
                                 init_vals=NULL,
                                 hyperparams=NULL,
                                 prev_covar=NULL,
                                 basis_functions=FALSE) {

  #check Y is a matrix 
  nn <- nrow(Y)
  categories <- apply(Y, 1, which.max)

  p <- ncol(X)
  K <- ncol(Y)
  r <- ncol(Psi)
  n_tp <- max(timepoints)

  if(!is.null(prev_covar)) {
    Z <- model.matrix(~ prev_covar:factor(timepoints) -1)
    Z <- Z[ , -c(2:(K+1))]
    n_gammas <- (K-1) * ncol(Z) 
  } else {
    n_gammas <- K-1
  }

  I_r <- Diagonal(r)

  default_init_vals <- list(beta = rep(0, p),
                            eta = matrix(rnorm(r*n_tp), nrow=r, ncol=n_tp),
                            sigma2_eta = 1,
                            sigma2_eta_1 = 1,
                            gamma = rnorm(n_gammas, mean=rep(0, n_gammas), sd=1),
                            phi = .5)

  init_vals <- validate_list_arg(default_init_vals, init_vals)

  default_hyperparams <- list(sigma2_beta = 3,
                              sigma2_gamma = 3,
                              a=.1,
                              b=.1)
  
  hyperparams <- validate_list_arg(default_hyperparams, hyperparams)

  beta_out <- matrix(NA, n_iter, p)
  eta_out <- array(NA, dim=c(n_iter, r, n_tp))
  gamma_out <- matrix(NA, n_iter, n_gammas)
  phi_out <- sigma2_eta_out <- sigma2_eta_1_out <-rep(NA, n_iter)

  I_g_sigg <- Diagonal(n_gammas)/hyperparams$sigma2_gamma
  I_p_sigb <- Diagonal(p)/hyperparams$sigma2_beta

  eta <- init_vals$eta
  beta <- init_vals$beta
  gamma <- init_vals$gamma
  sigma2_eta <- init_vals$sigma2_eta

  omega <- matrix(0, nn, K-1) #each column is the diagonal of Omega_k

  keep <- 1:nn
  G_ind <- matrix(0, nrow=nn, ncol=(K-1))
  G_ind[,1] <- 1
  kappa <- as.matrix(weights * (Y[,-K] - 0.5)) #column for each category
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
  timepoints <- timepoints[keep]
  if(!is.null(prev_covar)){
    Z <- Z[keep,]
    G_ind <- t(KhatriRao(t(G_ind), t(Z)))
  }

  gamma_ind <- G_ind %*% gamma
  X_beta <- X %*% beta

  sigma2_eta_1_shape <- hyperparams$a + r/2
  sigma2_eta_shape <- hyperparams$a + r*(n_tp-1)/2

  pb <- txtProgressBar(min=0, max=n_iter, style=3)
  for(i in 1:n_iter){
    Psi_eta <- (Psi %*% eta)[cbind(seq_along(timepoints), timepoints)]
    linpred <- as.numeric(gamma_ind - X_beta - Psi_eta)

    omega <- rpg.gamma(num=length(keep), h=weights, z=linpred)
    Omega <- Diagonal(n=length(keep), x=omega)
    kappa_omega <- kappa/omega

    denom <- sum(eta[, 1:(n_tp-1)] * eta[, 1:(n_tp-1)]) 
    mu_phi <- sum(eta[, 2:n_tp] * eta[, 1:(n_tp-1)])/denom
    sigma_phi <- sqrt(sigma2_eta/denom)
    phi <- phi_out[i] <- rtnorm(n=1, mean=mu_phi, sd=sigma_phi, a=-1, b=1)

    sigma2_eta_1 <- sigma2_eta_1_out[i] <- rinvgamma(n=1,
                                                     a=sigma2_eta_1_shape,
                                                     b=hyperparams$b + t(eta[,1])%*%eta[,1]/2)

    eta_diff <- eta[, 2:n_tp] - phi * eta[, 1:(n_tp-1)]
    sigma2_eta_scale <- hyperparams$b + sum(eta_diff * eta_diff)/2

    sigma2_eta <- sigma2_eta_1_out[i] <- rinvgamma(n=1,
                                                   a=sigma2_eta_shape,
                                                   b=sigma2_eta_scale)

    tXOmega <-t(X) %*% Omega 
    mu_beta <- tXOmega %*% (gamma_ind - kappa_omega - Psi_eta)
    prec_beta <- tXOmega %*% X + I_p_sigb
 
    tGOmega <- t(G_ind) %*% Omega
    mu_gamma <- tGOmega %*% (kappa_omega + X_beta + Psi_eta)
    prec_gamma <- tGOmega %*% G_ind + I_g_sigg

    mu_block <- c(as.vector(mu_beta), as.vector(mu_gamma))
    prec_block <- bdiag(prec_beta, prec_gamma)
    block <- rmvnorm_prec(n=1,
                          mean=solve(prec_block, mu_block),
                          prec=prec_block)
      
    beta <- beta_out[i, ] <- block[c(1:p)]
    gamma <- gamma_out[i, ] <- block[-c(1:p)]
    X_beta <- X %*% beta
    gamma_ind <- G_ind %*% gamma

    ids_1 <- which(timepoints == 1)
    Psi_1 <- Psi[ids_1, ]
    X_1 <- X[ids_1, ]
    Omega_1 <- Omega[ids_1, ids_1]
    kappa_omega_1 <- kappa_omega[ids_1]
    gamma_ind_1 <- gamma_ind[ids_1]
  
    tPsiOmega_1 <- t(Psi_1)%*%Omega_1
    prec_eta_1 <- tPsiOmega_1 %*%Psi_1 + (1/sigma2_eta_1 + phi^2/sigma2_eta)*I_r
    mu_eta_1 <- tPsiOmega_1%*%(gamma_ind_1 - kappa_omega_1 - X_1%*%beta) +
         phi/sigma2_eta*eta[,2]
    mu_eta_1 <- solve(prec_eta_1, mu_eta_1)
    eta[, 1] <- rmvnorm_prec(n=1,
                             mean = mu_eta_1,
                             prec=prec_eta_1)
    #sum to zero constraint
    #eta[, 1] <- eta[, 1] - mean(eta[, 1])

    for(tt in 2:(n_tp-1)){
      ids_t <- which(timepoints == tt)
      Psi_t <- Psi[ids_t, ]
      Omega_t <- Omega[ids_t, ids_t]
      X_t <- X[ids_t, ]
      kappa_omega_t <- kappa_omega[ids_t]
      gamma_ind_t <- gamma_ind[ids_t]
 
      tPsiOmega_t <- t(Psi_t)%*%Omega_t
      prec_eta_t <- tPsiOmega_t%*%Psi_t + (1+phi^2)/sigma2_eta * I_r
      mu_eta_t <- tPsiOmega_t%*%(gamma_ind_t - kappa_omega_t - X_t%*%beta) +
          phi/sigma2_eta*(eta[,(tt-1)] + eta[,(tt+1)])
      mu_eta_t <- solve(prec_eta_t, mu_eta_t)
      eta[, tt] <- rmvnorm_prec(n=1,
                                mean = mu_eta_t,
                                prec=prec_eta_t)
    }

    ids_T <- which(timepoints == n_tp)
    Psi_T <- Psi[ids_T,]
    Omega_T <- Omega[ids_T, ids_T]
    X_T <- X[ids_T,]
    kappa_omega_T <- kappa_omega[ids_T]
    gamma_ind_T <- gamma_ind[ids_T]
     
    tPsiOmega_T <- t(Psi_T)%*%Omega_T
    prec_eta_T <- tPsiOmega_T%*%Psi_T + I_r/sigma2_eta
    mu_eta_T <- tPsiOmega_T%*%(gamma_ind_T - kappa_omega_T - X_T%*%beta) +
         phi/sigma2_eta*eta[,(n_tp-1)]
    mu_eta_T <- solve(prec_eta_T, mu_eta_T)
    eta[, n_tp] <- rmvnorm_prec(n=1,
                                mean = mu_eta_T,
                                prec=prec_eta_T)

    eta_out[i,,] <- eta
      
    setTxtProgressBar(pb, i)
  }
  close(pb)
    
  return(list(posteriors =  list(eta=eta_out[-c(1:n_burn),, ],
                                 beta=beta_out[-c(1:n_burn), ],
                                 phi=phi_out[-c(1:n_burn)],
                                 sigma2_eta=sigma2_eta_out[-c(1:n_burn)],
                                 sigma2_eta_1=sigma2_eta_1_out[-c(1:n_burn)],
                                 gamma=gamma_out[-c(1:n_burn), ]),
              hyperparams=hyperparams))
}

