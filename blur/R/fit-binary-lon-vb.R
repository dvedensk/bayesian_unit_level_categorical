
#fit binary longitudinal with VB
fit_bin_lon_vb <- function(sigma2_beta=10000, Psi, X, Y, weights, # true_prec,
                         weeks, areas, epsilon, basis_funcs=NULL){#, true_phi){ 
  r <- ncol(Psi)
  n_weeks <- max(weeks)
  n_rt <- r*n_weeks
  p <- ncol(X)
  n_param <- n_rt + p

  I_p <- Diagonal(p)
  zero_p <- rep(0, p)
  I_r <- Diagonal(r)
  zero_r <- rep(0, r)

  A <- model.matrix(~ 0 + as.factor(weeks))
  B <- Psi# basis_funcs[areas,]#Psi#[areas,]
  PsiTime <- t(KhatriRao(t(A), t(B)))
  #indices for first and last timepoints in Eta
  #need these for sigma_eta and phi
  idx <- c(1 + n_weeks*(0:(r-1)),
           (1:r)*n_weeks)  

    ###
    ### Set initial values
    ###
    E_Beta <- rep(1, p) 
    X_Beta <- X %*% E_Beta
    sigma2_eta <- 1
    E_Eta <- Matrix(1, nrow=r, ncol=n_weeks)
    Sigma_Eta <- Diagonal(n=r*n_weeks) #block diag of all Sigma etas
    Sigma_Beta <- I_p
    a <- b <- 1
    E_sigma2_eta_inv <- E_sigma2_eta_1_inv <- 1
    a_tilde_1 <- a + r/2
    a_tilde <- a + r*(n_weeks-1)/2
    gam_tilde_1 <- gamma(a_tilde_1 + 1/2)/gamma(a_tilde_1)
    gam_tilde <- gamma(a_tilde + 1/2)/gamma(a_tilde)
    E_phi <- .5
    E_phi2 <- E_phi^2

    C <- Matrix(cbind(X, PsiTime))
    E_Zeta <- rep(1, p + n_rt)
    Sigma_Zeta <- bdiag(Sigma_Beta, Sigma_Eta)
    ###
    ### Start sampling
    ###
    kappa <- weights*(Y - .5)
    Psi_Eta <- rep(0, nrow(X_Beta))

    checkOld <- Inf
    checkVec <- c()
    iter <- 1 

    repeat{
     Psi_Eta <- PsiTime%*%as.vector(E_Eta)
        if(iter %% 100 == 0){print(paste0("iter = ", iter))}
        
     xi <- colSums(t(X)*(Sigma_Beta%*%t(X))) + 
	     colSums(t(PsiTime)*(Sigma_Eta%*%t(PsiTime))) + 
	     (X_Beta + Psi_Eta)^2
     xi <- as.numeric(sqrt(xi))

     omega <- (weights*0.5/xi)*tanh(0.5*xi)
     E_Omega <- Diagonal(x=omega)
     kappa_omega <- kappa/omega
       
     #phi
     phi_lower <- -1
     phi_upper <- 1

     numer <- sum(E_Eta[,-1]*E_Eta[,-n_weeks])
     denom <- sum(E_Eta[,-n_weeks] * E_Eta[,-n_weeks]) +
                 sum(diag(Sigma_Eta[1:(r*(n_weeks-1)), 1:(r*(n_weeks-1))]))
     mu_phi <- numer/denom
     sigma2_phi <- 1/(E_sigma2_eta_inv * denom)
     sigma_phi <- sqrt(sigma2_phi)

     alpha_phi <- (phi_lower - mu_phi)/sigma_phi
     beta_phi <- (phi_upper - mu_phi)/sigma_phi
     pnorm_diff <- pnorm(beta_phi) - pnorm(alpha_phi)
     dnorm_diff <- dnorm(beta_phi) - dnorm(alpha_phi)
     E_phi <- mu_phi - sigma_phi * dnorm_diff/pnorm_diff
     Var_phi <- sigma2_phi *
         (1 - (beta_phi*dnorm(beta_phi) - alpha_phi*dnorm(alpha_phi))/pnorm_diff -
          (dnorm_diff/pnorm_diff)^2)
     E_phi2 <- 1 - (beta_phi*dnorm(beta_phi) - alpha_phi*dnorm(alpha_phi))/pnorm_diff

     E_phi <- mu_phi
     Var_phi <- sigma2_phi
     E_phi2 <- Var_phi + E_phi^2
#     E_phi <- true_phi
#     E_phi2 <- true_phi^2

     ##sample sigma_eta_1
     b_tilde_1 <- b + 0.5*(sum(E_Eta[,1]*E_Eta[,1]) + sum(diag(Sigma_Eta[1:r,1:r])))
     E_sigma2_eta_1_inv <- a_tilde_1/b_tilde_1

     #sample sigma^2_eta
     tmp1 <- sum(E_Eta[,-1] * E_Eta[,-1]) + sum(diag(Sigma_Eta[-(1:r), -(1:r)]))
     tmp2 <- -2*E_phi * sum(E_Eta[,-1]*E_Eta[,-n_weeks])
     tmp3 <- E_phi2 * (sum(E_Eta[,-n_weeks] * E_Eta[,-n_weeks]) +
     		  sum(diag(Sigma_Eta[1:(r*(n_weeks-1)), 1:(r*(n_weeks-1))])))
     b_tilde <- 0.5*(tmp1+tmp2+tmp3)
     E_sigma2_eta_inv <- a_tilde/b_tilde
#     E_sigma2_eta_inv <- true_prec

     #beta
     tXOmega <-t(X) %*% E_Omega
     Prec_Beta <- tXOmega %*% X + 1/sigma2_beta * I_p
     Sigma_Beta <- solve(Prec_Beta)
     E_Beta <- Sigma_Beta%*%tXOmega %*% (kappa_omega - Psi_Eta) 
     X_Beta <- X %*% E_Beta         

     #etas
     ids_1 <- which(weeks==1)
     Psi_1 <- Psi[ids_1,]
     X_1 <- X[ids_1,]
     E_Omega_1 <- E_Omega[ids_1, ids_1]
     kappa_omega_1 <- kappa_omega[ids_1]

     tPsiOmega_1 <- t(Psi_1) %*% E_Omega_1
     Prec_Eta <- tPsiOmega_1 %*% Psi_1 + (E_sigma2_eta_1_inv + E_phi2*E_sigma2_eta_inv)*I_r
     Sigma_Eta[1:r,1:r] <- solve(Prec_Eta)
     E_Eta[,1] <- as.vector(Sigma_Eta[1:r,1:r]%*%(tPsiOmega_1 %*% (kappa_omega_1 - X_1%*%E_Beta) + E_phi * E_sigma2_eta_1_inv * E_Eta[,2]))
     E_Eta[,1] <- E_Eta[,1]-mean(E_Eta[,1])

     ##sample eta[,2:(n_weeks-1)]
     for(tt in 2:(n_weeks-1)){
       ids_t <- which(weeks==tt)
       Psi_t <- Psi[ids_t,]
       E_Omega_t <- E_Omega[ids_t, ids_t]
       X_t <- X[ids_t,]
       kappa_omega_t <- kappa_omega[ids_t]
       id_r <- (1:r)+(tt-1)*r

       tPsiOmega_t <- t(Psi_t) %*% E_Omega_t
       Prec_Eta_t <- tPsiOmega_t%*%Psi_t + (1+E_phi2)*E_sigma2_eta_inv * I_r
       Prec_Eta <- bdiag(Prec_Eta, Prec_Eta_t)
       Sigma_Eta[id_r,id_r] <-solve(Prec_Eta_t)
       E_Eta[,tt] <- as.vector( Sigma_Eta[id_r, id_r]%*%
                     (t(Psi_t)%*%E_Omega_t%*%(kappa_omega_t - X_t%*%E_Beta) +
                          E_phi*E_sigma2_eta_inv*(E_Eta[,(tt-1)] + E_Eta[,(tt+1)])))
     }

     #sample eta[,n_weeks]
     ids_T <- which(weeks==n_weeks)
     Psi_T <- Psi[ids_T,]
     E_Omega_T <- E_Omega[ids_T, ids_T]
     X_T <- X[ids_T,]
     kappa_omega_T <- kappa_omega[ids_T]

     id_r <- (1:r)+r*(n_weeks-1)
     tPsiOmega_T <- t(Psi_T) %*% E_Omega_T
     Prec_Eta_T <- tPsiOmega_T %*% Psi_T + I_r*E_sigma2_eta_inv
     Prec_Eta <- bdiag(Prec_Eta, Prec_Eta_T)
     Sigma_Eta[id_r, id_r] <- solve(Prec_Eta_T)
     E_Eta[,n_weeks] <- Sigma_Eta[id_r,id_r]%*%
                             (t(Psi_T)%*%E_Omega_T%*%(kappa_omega_T - X_T%*%E_Beta) +
                              E_phi*E_sigma2_eta_inv * E_Eta[,(n_weeks-1)])

#     Psi_Eta <- (Psi %*% E_Eta)[cbind(seq_along(weeks), weeks)]

     PPrec <- bdiag(Prec_Beta, Prec_Eta)
     Sigma_Zeta <- bdiag(Sigma_Beta, Sigma_Eta)
     E_Zeta <- c(as.vector(E_Beta), as.vector(E_Eta))
        
     checkNew  <- 0.5*(n_param) +
         0.5*determinant(Sigma_Zeta, logarithm = T)$modulus  -
         0.5*t(E_Zeta)%*%PPrec%*%(E_Zeta) + 
 #       sum(weights*(Y-0.5)*as.numeric(D%*%MUbu) + 
         sum(kappa*as.numeric(C%*%E_Zeta) + 
         weights*log(plogis(xi)) - 0.5*weights*xi) - 
         0.5*sum(diag(PPrec %*% E_Zeta)) - log(b_tilde)
     checkVec <- c(checkVec, as.numeric(checkNew))
     
     if (as.logical(abs(checkOld - checkNew) < epsilon)) break
     if (iter > 1000) break
     iter <- iter + 1
     checkOld <- checkNew
  }
  return(list(E_Zeta=E_Zeta, Sigma_Zeta=Sigma_Zeta, E_phi=E_phi, E_sigma2_eta_1_inv=E_sigma2_eta_1_inv,
              E_sigma2_eta_inv=E_sigma2_eta_inv, Elbo=checkVec))
}

