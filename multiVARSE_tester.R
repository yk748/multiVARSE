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
  # Compute critical values and make deicisions:
  if (type == "nullity"|type == "homogeneity"){
    # ------------------------------------------- #
    # Compute D and M:
    if (type == "nullity"){
      D <- diag(1,K)
      df <- K
    }else if (type == "homogeneity"){
      D <- matrix(0,(K-1),K)
      for (k in 1:(K-1)){
        D[k,k] <- 1
        D[k,(k+1)] <- -1
      } 
      df <- K-1
    }
    M <- solve(diag(sqrt(TT),K))
    
    # ------------------------------------------- #
    # Compute test statistics:
    test_ij <- array(NA,dim=c(d,(d*p)))
    for (i in 1:d){
      for (j in 1:(d*p)){
        idx_single <- ((j-1)*d+i)
        S_ij <- D %*% M %*% diag(V_ij_mat[i,j,]) %*% M %*% t(D)
        Db_ij <- D %*% mapply(k=1:K,function(k)debias[[k]]$beta_tilde[idx_single])
        test_ij[i,j] <-  t(Db_ij) %*% solve(S_ij) %*% Db_ij
      }
    }
    critical <- qchisq((1-alpha),df)
    decision <- array(test_ij >= critical,dim=c(d,d,p))
    
  }else if (type == "significance"){
    
    idx_arr <- t(array(c(1:(d^2*p)),dim=c(d*p,d)))
    V_ij_tot <- N_ij_tot <- array(NA,dim=c(d,d*p))
    test_ij <- array(NA,dim=c(d,d*p))
    for (i in 1:d){
      for (j in 1:(d*p)){
        V_ij_sum <- 0
        subj_idx <- identified$ind_contribute[[idx_arr[i,j]]]
        for (k in 1:length(subj_idx)){
          V_ij_sum <- V_ij_sum + abs(V_ij_mat[i,j,k])
        }
        V_ij_tot[i,j] <- V_ij_sum/length(subj_idx)^2
        N_ij_tot[i,j] <- mean(TT[subj_idx])
      }
    }
    test_ij <- sqrt(N_ij_tot)*matrix(identified$alpha_0_hat,d,(d*p),byrow=FALSE)/sqrt(V_ij_tot)
    critical <- qnorm(1-alpha/2)
    decision <- array(abs(test_ij) >= critical,dim=c(d,d,p))
  }
  reject_idx <- which(decision==TRUE,arr.ind=TRUE)[,-3]
  nonreject_idx <- which(decision==FALSE,arr.ind=TRUE)[,-3]

  # ---------------------------------------------------- #
  # Return output:
  output <- list(sigma2_hat=sigma2_hat,Omega_hat=Omega_hat,
                 type = type,
                 test_stat = test_ij,
                 critical = critical,
                 decision = decision,
                 reject_idx = reject_idx,
                 nonreject_idx = nonreject_idx)
  return(output)
}