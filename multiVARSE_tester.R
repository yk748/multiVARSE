###################################################################################
# The following is the code for performing three hypothesis tests 
# and their multiple test procedure explained in Sections 4.2 and 4.3.
###################################################################################
multiVARSE_tester <- function(X_t,debias,opt,identified,type,alpha,K,d,p,TT){
  
  X_t_scaled <- list()
  for (k in 1:K){
    X_t_scaled[[k]] <- t(scale(t(X_t[[k]])))
  }
  
  # ------------------------------------------- #
  # Compute sigma_X^2_hat:
  sigma2_hat <- array(NA,dim=c(d,K))
  for (k in 1:K){
    YY <- t(X_t_scaled[[k]][,TT[k]:(p+1)])
    if ( p == 1 ){
      XX <- t(X_t_scaled[[k]][,(TT[k]-1):p])
    }else{
      tmp_X <- t(X_t_scaled[[k]][,(TT[k]-1):p])
      for (l in 2:p){
        tmp_X <- cbind(tmp_X,t(X_t_scaled[[k]][,(TT[k]-l):(p-l)]))
      }
      XX <- tmp_X
    }
    
    for (i in 1:d){
      idx <- ((i-1)*d+1):(i*d)
      sigma2_hat[i,k] <- sum((YY[,i] - XX%*%opt[[k]]$beta_hat[idx])^2)/TT[k]
    }
  }
  
  # ------------------------------------------- #
  # Compute Sigma_X_hat and Omega_X_hat:
  Cov_X <- list()
  N_k <- vector("numeric",K)
  Omega_hat <- list()
  for (k in 1:K){
    N_k[k] <- TT[k] - p
    Cov_X[[k]] <- X_t_scaled[[k]] %*% t(X_t_scaled[[k]])/TT[k]
    Omega_hat[[k]] <- debias[[k]]$Theta_hat %*% Cov_X[[k]] %*% t(debias[[k]]$Theta_hat)
  }
  
  # ------------------------------------------- #
  # Compute V_hat:
  V_ij_mat <- array(NA,dim=c(d,(d*p),K))
  for (i in 1:d){
    for (j in 1:(d*p)){
      V_ij_mat[i,j,] <- diag(diag(sigma2_hat[i,]) %*% diag(mapply(k=1:K,function(k)Omega_hat[[k]][j,j])))
    }
  }
  
  # ------------------------------------------- #
  # Compute critical values and make decisions:
  if (type == "homogeneity"){
    # ------------------------------------------- #
    # pre-steps:
    pairs <- combn(1:K, 2)
    n_sim <- 500
    
    # ------------------------------------------- #
    # Compute test statistics:
    test_ij <- array(NA,dim=c(d,(d*p)))
    for (i in 1:d){
      for (j in 1:(d*p)){
        idx_single <- ((j-1)*d+i)
        beta_tilde_ij <- mapply(k=1:K,function(k)debias[[k]]$beta_tilde[idx_single])
        
        result <- compute_M_lambda(V_ij_mat[i,j,], TT)
        lambda <- result$lambda
        M <- result$M
        Z <- matrix(rnorm(M * n_sim), nrow = M, ncol = n_sim)
        Q <- sqrt(colSums(lambda*Z^2))
        
        T_ij <- vector("numeric",dim(pairs)[2])
        for (l in 1:dim(pairs)[2]){
          k1 <- pairs[1,l]
          k2 <- pairs[2,l]
          T_ij[l] <- (beta_tilde_ij[k1] - beta_tilde_ij[k2])/sqrt(V_ij_mat[i,j,k1]/TT[k1] + V_ij_mat[i,j,k2]/TT[k2])
        }
        test_ij[i,j] <- qnorm(pmax(pmin(sum(norm(T_ij,"2") >= Q)/n_sim, 1-1e-15), 1e-15))
      }
    }
    
  }else if (type == "nullity"){
    # ------------------------------------------- #
    # Compute test statistics:
    test_ij <- array(NA,dim=c(d,(d*p)))
    for (i in 1:d){
      for (j in 1:(d*p)){
        idx_single <- ((j-1)*d+i)
        beta_tilde_ij <- mapply(k=1:K,function(k)debias[[k]]$beta_tilde[idx_single])
        chi_sq_ij <- sum( (sqrt(TT)*beta_tilde_ij/sqrt(V_ij_mat[i,j,]))^2 )
        test_ij[i,j] <- qnorm(pchisq(chi_sq_ij, df = K))
      }
    }
    
  }else if (type == "significance"){
    
    idx_arr <- t(array(c(1:(d^2*p)),dim=c(d*p,d)))
    V_ij_tot <- N_ij_tot <- array(NA,dim=c(d,d*p))
    test_ij <- array(NA,dim=c(d,d*p))
    for (i in 1:d){
      for (j in 1:(d*p)){
        V_ij_sum <- 0
        subj_idx <- identified$ind_contribute[[idx_arr[i,j]]]
        if (length(subj_idx)!=0){
          for (k in 1:length(subj_idx)){
            V_ij_sum <- V_ij_sum + abs(V_ij_mat[i,j,k])
          }
          V_ij_tot[i,j] <- V_ij_sum/length(subj_idx)^2
          N_ij_tot[i,j] <- mean(TT[subj_idx])
        }else{
          V_ij_tot[i,j] <- NA
          N_ij_tot[i,j] <- NA
        }
      }
    }
    test_ij <- sqrt(N_ij_tot)*matrix(identified$alpha_0_hat,d,(d*p),byrow=FALSE)/sqrt(V_ij_tot)
  }
  
  # ------------------------------------------- #
  # Fixed critical value (single test, no FDR control):
  critical <- qnorm(1 - alpha/2)
  decision <- array(abs(test_ij) >= critical, dim = c(d,d,p))
  reject_idx <- which(decision == TRUE,  arr.ind = TRUE)[,-3]
  nonreject_idx <- which(decision == FALSE, arr.ind = TRUE)[,-3]
  
  # ------------------------------------------- #
  # BH-type FDR threshold (eqref{e:thres_sig}):
  # hat_zeta = inf { 0 <= zeta <= sqrt(2 log |S|) :
  #   2 * Phibar(zeta) * |S| / max(1, #{|Z| >= zeta}) <= alpha }
  Z_vec <- as.vector(test_ij)          # all Z_{(i,j)}, |S| = d*(d*p)
  all_idx <- which(!is.na(matrix(Z_vec, d, d*p)), arr.ind = TRUE)
  S_size <- length(Z_vec)
  zeta_upper <- sqrt(2*log(S_size))
  
  # Search over observed |Z| values in decreasing order
  # (denominator only changes at observed values)
  zeta_candidates <- sort(abs(Z_vec[!is.na(Z_vec)]), decreasing = TRUE)
  zeta_candidates <- zeta_candidates[zeta_candidates <= zeta_upper]
  
  hat_zeta <- zeta_upper  # default: no interior solution found
  for (zeta in zeta_candidates) {
    num_reject <- sum(abs(Z_vec) >= zeta, na.rm = TRUE)
    ratio <- 2*pnorm(zeta, lower.tail=FALSE)*S_size/max(1, num_reject)
    if (ratio <= alpha){
      hat_zeta <- zeta
      break  # take the largest zeta satisfying condition = infimum
    }
  }
  
  # Controlled decisions:
  decision_controlled <- array(abs(test_ij) >= hat_zeta, dim = c(d, d, p))
  reject_idx_controlled    <- which(decision_controlled == TRUE,  arr.ind = TRUE)[,-3]
  nonreject_idx_controlled <- which(decision_controlled == FALSE, arr.ind = TRUE)[,-3]
  
  # ---------------------------------------------------- #
  # Return output:
  output <- list(
    sigma2_hat    = sigma2_hat,
    Omega_hat     = Omega_hat,
    type          = type,
    test_stat     = test_ij,       # Z_mat for evaluator
    all_idx       = all_idx,       # all (i,j) in S
    # Single test (fixed critical value):
    critical      = critical,
    decision      = decision,
    reject_idx    = reject_idx,
    nonreject_idx = nonreject_idx,
    # Multiple test (BH-type FDR control):
    hat_zeta                 = hat_zeta,
    decision_controlled      = decision_controlled,
    reject_idx_controlled    = reject_idx_controlled,
    nonreject_idx_controlled = nonreject_idx_controlled
  )
  return(output)
}


#################################################################################
# Find M and lambda:
compute_M_lambda <- function(V_ij, TT, tol = 1e-10) {
  
  K <- length(TT)
  pairs <- which(lower.tri(matrix(0,K,K)), arr.ind=TRUE)
  pairs <- pairs[order(pairs[, 2], pairs[, 1]),]
  pairs <- t(apply(pairs, 1, function(x) sort(x)))
  pairs <- unique(pairs)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), ]
  n_pairs <- nrow(pairs)
  
  # ------------------------------------------- #  
  # Fill in Sigma_{(i,j)} using formula
  denom <- function(k, l) {
    sqrt(V_ij[k]/TT[k] + V_ij[l]/TT[l])
  }
  numer <- function(k,l,kp,lp) {
    (V_ij[k]/TT[k]*(k==kp) + V_ij[l]/TT[l]*(l==lp) 
     - V_ij[k]/TT[k]*(k==lp) - V_ij[l]/TT[l]*(l==kp))
  }
  
  Sigma_ij <- matrix(0, nrow=n_pairs,ncol=n_pairs)
  for (a in 1:n_pairs) {
    for (b in 1:n_pairs) {
      k  <- pairs[a,1]; l  <- pairs[a,2]; kp <- pairs[b,1]; lp <- pairs[b,2]
      Sigma_ij[a, b] <- numer(k,l,kp,lp) / (denom(k,l)*denom(kp,lp))
    }
  }
  
  # ------------------------------------------- #
  # Eigen decomposition of Sigma_{(i,j)} and get M and lambda_1,...,lambda_M
  eig    <- eigen(Sigma_ij, symmetric = TRUE)
  lambda <- eig$values
  pos_idx <- lambda > tol
  lambda_pos <- lambda[pos_idx]
  M <- sum(pos_idx)
  
  output <- list( M = M,                  
                  lambda = lambda_pos,        
                  Sigma_ij = Sigma_ij,         
                  pairs = pairs,           
                  lambda_all = lambda)
  return(output)
}
