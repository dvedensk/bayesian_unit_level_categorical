#' This function fits the longitudinal ordinal model via variational Bayes
#' @import Matrix
#' @param Y A vector of binary survey responses
#' @param X The input covariate matrix, one row per sample unit
#' @param Psi The matrix of spatial basis functions, one row per sample unit
#' @param timepoints A vector of the time at which measurements were taken
#' @param weights The vector of survey weights
#' @param epsilon is the convergence tolerance
#' @param n_samples number of samples to draw from variational posterior
#' @param init_vals Initial values
#' @param hyperparams Fixed hyperparameters for fixed effects variance and inverse gamma priors
#' @param prev_covar List of previous responses if longitudinal effect is desired

fit_ordinal_lon_vb <- function(Y,
                               X,
                               Psi,
                               weights,
                               timepoints,
                               epsilon=0.001,
                               n_samples,
                               init_vals=NULL,
                               hyperparams=NULL,
                               prev_covar=NULL) {

  nn <- nrow(Y)
  categories <- apply(Y, 1, which.max)
  areas <- apply(Psi, 1, which.max)

  p <- ncol(X)
  K <- ncol(Y)
  r <- ncol(Psi)
  n_tp <- max(timepoints)

  if(!is.null(prev_covar)) {
    Z <- model.matrix(~ prev_covar:factor(timepoints) -1)
    Z <- Z[,-c(2:(K+1))]
    n_gammas <- (K-1) * ncol(Z) 
  } else{
    n_gammas <- K-1
  }

  I_r <- Diagonal(r)
  I_p <- Diagonal(p)
  I_g <- Diagonal(n_gammas)
  
  MLEs <- MASS::polr(ordered(categories) ~ X)
    
  default_init_vals <- list(beta = MLEs$coefficients,
#                            eta = matrix(rnorm(r*n_tp), nrow=r, ncol=n_tp),
                            eta = matrix(1, nrow=r, ncol=n_tp),
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
    
  I_g_sigg <- I_g/hyperparams$sigma2_gamma
  I_p_sigb <- I_p/hyperparams$sigma2_beta
  
  E_Beta <- init_vals$beta
  E_Gamma <- init_vals$gamma
  omega <- matrix(1, nn, K-1) #each column is the diagonal of Omega_k
  E_sigma2_eta_1_inv <- 1/init_vals$sigma2_eta_1
  E_sigma2_eta_inv <- 1/init_vals$sigma2_eta
  E_Eta <- init_vals$eta
  Sigma_Eta <- Diagonal(n=r*n_tp) #block diag of all Sigma etas
  Sigma_Beta <- I_p
  Sigma_Gamma <- I_g

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
  Psi <- Psi[keep,]
  weights <- weights[keep]
  timepoints <- timepoints[keep]
  areas <- areas[keep]
  if(!is.null(prev_covar)){
    Z <- Z[keep,]
    G_ind <- t(KhatriRao(t(G_ind), t(Z)))
  }

  A <- model.matrix(~ 0 + as.factor(timepoints))
  B <- Psi #basis_funcs[areas,]
  PsiTime <- t(KhatriRao(t(A), t(B)))

  a_tilde_1 <- hyperparams$a + r/2
  a_tilde <- hyperparams$a + r*(n_tp-1)/2
  
  checkOld <- Inf

  X_Beta <- X%*%E_Beta
  gamma_ind <- G_ind %*% E_Gamma
  check_new  <- check_diff <- Inf

  Psi_Eta <- PsiTime%*%as.vector(E_Eta)

  repeat{
   xi <- colSums(t(G_ind)*(Sigma_Gamma%*%t(G_ind))) +
           colSums(t(X)*(Sigma_Beta%*%t(X))) + #tmpPsiEta + 
           colSums(t(PsiTime)*(Sigma_Eta%*%t(PsiTime))) +
           (gamma_ind - X_Beta - Psi_Eta)^2
   xi <- as.numeric(sqrt(xi))

   omega <- (weights*0.5/xi)*tanh(0.5*xi)
   E_Omega <- Diagonal(x=omega)
   kappa_omega <- kappa/omega
  
   ##Phi
   phi_lower <- -1
   phi_upper <- 1

   numer <- sum(E_Eta[,-1]*E_Eta[,-n_tp])
   denom <- sum(E_Eta[,-n_tp] * E_Eta[,-n_tp]) +
                sum(diag(Sigma_Eta[1:(r*(n_tp-1)), 1:(r*(n_tp-1))]))
   mu_phi <- numer/denom
   sigma2_phi <- 1/(E_sigma2_eta_inv * denom)
   sigma_phi <- sqrt(sigma2_phi)

   alpha <- (phi_lower - mu_phi)/sigma_phi
   beta <- (phi_upper - mu_phi)/sigma_phi
   pnorm_diff <- pnorm(beta) - pnorm(alpha)
   dnorm_diff <- dnorm(beta) - dnorm(alpha)
   E_phi <- mu_phi - sigma_phi * dnorm_diff/pnorm_diff 
   Var_phi <- sigma2_phi *
        (1 - (beta*dnorm(beta) - alpha*dnorm(alpha))/pnorm_diff -
          (dnorm_diff/pnorm_diff)^2)

   E_phi2 <- Var_phi + E_phi^2

   ##sigma eta 1 and sigma eta
   b_tilde_1 <- hyperparams$b + 0.5*(sum(E_Eta[,1]*E_Eta[,1]) + sum(diag(Sigma_Eta[1:r,1:r])))
   E_sigma2_eta_1_inv <- a_tilde_1/b_tilde_1
      
   #sample sigma^2_eta
   #trace for all but time 1:
   tmp1 <-  sum(E_Eta[,-1] * E_Eta[,-1]) + sum(diag(Sigma_Eta[-(1:r), -(1:r)]))
   tmp2 <- -2*E_phi * sum(E_Eta[,-1]*E_Eta[,-n_tp])
   #trace for all but time T
   tmp3 <-   E_phi2 * (sum(E_Eta[,-n_tp] * E_Eta[,-n_tp]) +
                           sum(diag(Sigma_Eta[1:(r*(n_tp-1)), 1:(r*(n_tp-1))])))
   b_tilde <- hyperparams$b + 0.5*(tmp1 + tmp2 + tmp3)
   E_sigma2_eta_inv <- a_tilde/b_tilde
     
   ##Beta
   tXOmega <-t(X) %*% E_Omega 
   Prec_Beta <- tXOmega %*% X + I_p_sigb
   Sigma_Beta <- solve(Prec_Beta)
   E_Beta <- Sigma_Beta%*%tXOmega %*% (gamma_ind - kappa_omega - Psi_Eta)
   X_Beta <- X %*% E_Beta

   ##Gamma
   tGOmega <- t(G_ind) %*% E_Omega
   Prec_Gamma <- tGOmega %*% G_ind + I_g_sigg
   Sigma_Gamma <- solve(Prec_Gamma)
   E_Gamma <- Sigma_Gamma%*% tGOmega %*% (kappa_omega + X_Beta + Psi_Eta)
   gamma_ind <- G_ind %*% E_Gamma

   ##Etas
   ids_1 <- which(timepoints==1)
   Psi_1 <- Psi[ids_1,]
   X_1 <- X[ids_1,]
   E_Omega_1 <- E_Omega[ids_1, ids_1]
   kappa_omega_1 <- kappa_omega[ids_1]
   gamma_ind_1 <- gamma_ind[ids_1]
  
   tPsiOmega_1 <- t(Psi_1) %*% E_Omega_1
   Prec_Eta <- tPsiOmega_1 %*% Psi_1 + (E_sigma2_eta_1_inv + E_phi2*E_sigma2_eta_inv)*I_r
   Sigma_Eta[1:r,1:r] <- solve(Prec_Eta)
   E_Eta[,1] <- as.vector(Sigma_Eta[1:r,1:r]%*%(tPsiOmega_1 %*% (gamma_ind_1 -
                                                                 kappa_omega_1 -
                                                                 X_1%*%E_Beta) +
                                                E_phi * E_sigma2_eta_1_inv * E_Eta[,2]))
#   E_Eta[,1] <- E_Eta[,1] - mean(E_Eta[,1])

   for(tt in 2:(n_tp-1)){
     ids_t <- which(timepoints==tt)
     Psi_t <- Psi[ids_t,]
     E_Omega_t <- E_Omega[ids_t, ids_t]
     X_t <- X[ids_t,]
     kappa_omega_t <- kappa_omega[ids_t]
     gamma_ind_t <- gamma_ind[ids_t]
     id_r <- (1:r)+(tt-1)*r

     tPsiOmega_t <- t(Psi_t) %*% E_Omega_t
     Prec_Eta_t <- tPsiOmega_t%*%Psi_t + (1+E_phi2)*E_sigma2_eta_inv * I_r
     Prec_Eta <- bdiag(Prec_Eta, Prec_Eta_t)
     Sigma_Eta[id_r,id_r] <-solve(Prec_Eta_t)
     E_Eta[,tt] <-as.vector( Sigma_Eta[id_r, id_r] %*% (t(Psi_t)%*%E_Omega_t%*%(gamma_ind_t -
                                                                                kappa_omega_t -
                                                                                X_t%*%E_Beta) +
                                                        E_phi*E_sigma2_eta_inv*(E_Eta[,(tt-1)] +
                                                                                E_Eta[,(tt+1)])))
     }

     ids_T <- which(timepoints==n_tp)
     Psi_T <- Psi[ids_T,]
     E_Omega_T <- E_Omega[ids_T, ids_T]
     X_T <- X[ids_T,]
     kappa_omega_T <- kappa_omega[ids_T]
     gamma_ind_T <- gamma_ind[ids_T]

     id_r <- (1:r)+r*(n_tp-1)
     tPsiOmega_T <- t(Psi_T) %*% E_Omega_T
     Prec_Eta_T <- tPsiOmega_T %*% Psi_T + I_r*E_sigma2_eta_inv
     Prec_Eta <- bdiag(Prec_Eta, Prec_Eta_T)
     Sigma_Eta[id_r, id_r] <- solve(Prec_Eta_T)
     E_Eta[, n_tp] <- as.vector(Sigma_Eta[id_r,id_r] %*% as.vector((t(Psi_T)%*%E_Omega_T%*%(gamma_ind_T -
                                                                                         kappa_omega_T -
                                                                                         X_T%*%E_Beta) +
                                                                 E_phi*E_sigma2_eta_inv * E_Eta[,(n_tp-1)])))

     Psi_Eta <- PsiTime%*%as.vector(E_Eta)

     PPrec <- bdiag(Prec_Gamma, Prec_Beta, Prec_Eta)
     Sigma_Zeta <- bdiag(Sigma_Gamma, Sigma_Beta, Sigma_Eta)
     E_Zeta <-c(as.vector(E_Gamma), as.vector(E_Beta), as.vector(E_Eta))

     MUbu <-  Sigma_Beta%*%t(X)%*%kappa
     checkNew  <- 0.5*(p + r*n_tp + n_gammas) +
          0.5*determinant(Sigma_Zeta, logarithm = T)$modulus 

    if (as.logical(abs(checkOld - checkNew) < epsilon)) break
    checkOld <- checkNew
  }
  variational_samples <- rmvnorm(n=n_samples,
                                 mean=E_Zeta,
                                 covar=Sigma_Zeta)

  gammas <- t(variational_samples[1:n_gammas, ])
  betas <- t(variational_samples[(n_gammas+1):(n_gammas+p), ])
  etas <- t(variational_samples[-c(1:(n_gammas+p)), ])
  etas <- array(etas, dim=c(n_samples, r, n_tp))

  sigma2_eta_1s <- rinvgamma(n=n_samples, alpha=a_tilde_1, beta=b_tilde_1)
  sigma2_etas <- rinvgamma(n=n_samples, alpha=a_tilde, beta=b_tilde)
  phis <- rtnorm(n=n_samples, a=-1, b=1, mean=E_phi, sd=sqrt(Var_phi))
 
  return(list(posteriors = list(test_mu=E_Zeta,
                                E_phi=E_phi,
                                beta=betas,
                                eta=etas,
                                gamma=gammas,
                                phi=phis,
                                sigma2_eta_1=sigma2_eta_1s,
                                sigma2_eta=sigma2_etas),
              hyperparams=hyperparams))    
}
