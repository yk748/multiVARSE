multiVARSE_cv <- function(X_t,debias,opt,K,d,p,TT,thres_type,nfold=5,
                          grid_eta=10,grid_thres_0=10,grid_thres_k=5){
  
  X_t_scaled <- list()
  kappa <- rep(NA,K)
  for (k in 1:K){
    X_t_scaled[[k]] <- t(scale(t(X_t[[k]])))
    
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
    vec_YY <- as.vector(YY)
    ZZ <- bdiag(replicate(d,XX,simplify=FALSE))
    E <- matrix(vec_YY - ZZ %*% opt[[k]]$beta_hat,ncol=d,byrow=FALSE)
    kappa[k] <- max(diag(t(E) %*% E / TT[k]))/min(diag(t(E) %*% E / TT[k]))
  }
  
  # ------------------------------------------- #
  # Setup:
  eta_min <- max(mapply(k=1:K,function(k)min(abs(debias[[k]]$beta_tilde))))
  eta_max <- min(mapply(k=1:K,function(k)max(abs(debias[[k]]$beta_tilde))))
  eta_grid <- seq(eta_min,eta_max,length=grid_eta)
  
  e_0 <- seq(0.1,1,length=grid_thres_0)
  c_k <- seq(0.5,1,length=grid_thres_k)
  thres_0_grid <- max(kappa)*sqrt(log(d^2*p)/(e_0*K*min(TT)))
  thres_k_grid <- mapply(m=1:grid_thres_k,function(m)c_k[m]*kappa*sqrt(log(d^2*p)/(TT))) # K x grid_thres_k
  
  loss_eta_thres <- array(NA,dim=c(grid_eta,grid_thres_0,grid_thres_k))
  
  # ------------------------------------------- #
  # Main loop:
  list_block_error <- list()
  for (l in 1:nfold){
    # ------------------------------------------- #
    # Training and test data for fixed fold 
    debiased_train <- list(); 
    TT_train_vec <- vector("numeric",K)
    vec_YY_test <- ZZ_test <- list()
    for (k in 1:K){
      # ------------------------------------------- #
      # Indexing:
      block_size <- round(TT[k]/nfold)
      if (l != nfold){
        idx_test <- (block_size*(l-1)+1):(block_size*l)
      }else{
        idx_test <- (block_size*(l-1)+1):TT[k]
      }
      idx_train <- c(1:TT[k])[-idx_test]
      
      # ------------------------------------------- #
      # Estimating & debiasing with train data:
      X_t_train <- X_t_scaled[[k]][,idx_train]
      TT_train_vec[k] <- dim(X_t_train)[2]
      opt <- multiVARSE_optimizer(X_t_train,d,p,TT_train_vec[k])
      debiased_train[[k]] <- multiVARSE_debiasor(opt$beta_hat,
                                                 X_t_train,d,p,TT_train_vec[k])
      
      # ------------------------------------------- #
      # Get YY and XX from test data:
      X_t_test <- X_t_scaled[[k]][,idx_test]
      TT_test <- dim(X_t_test)[2]
      YY_test <- t(X_t_test[,TT_test:(p+1)])
      if ( p == 1 ){
        XX_test <- t(X_t_test[,(TT_test-1):p])
      }else{
        tmp_X <- t(X_t_test[,(TT_test-1):p])
        for (l in 2:p){
          tmp_X <- cbind(tmp_X,t(X_t_test[,(TT_test-l):(p-l)]))
        }
        XX_test <- tmp_X
      }
      vec_YY_test[[k]] <- as.vector(YY_test)
      ZZ_test[[k]] <- bdiag(replicate(d,XX_test,simplify=FALSE))
      
      cat(k,"th subject done.\n")
    }
    
    # ------------------------------------------- #
    # Compute loss functions with fixed eta and thresholds:
    for (i in 1:grid_eta){
      for (j in 1:grid_thres_0){
        for (m in 1:grid_thres_k){
          identified_ijm <- multiVARSE_identifier(debiased_train,
                                                  eta_grid[i],thres_0_grid[j],thres_k_grid[,m],
                                                  thres_type=thres_type,
                                                  K,d,p,TT_train_vec)
          
          cat(i,"th eta,",j,"th thres0",m,"th thres k considered.\n")
          tmp_sum_ijm <- 0
          for (k in 1:K){
            TT_test <- length(vec_YY_test[[k]])
            tmp_sum_ijm <- tmp_sum_ijm + sum((vec_YY_test[[k]]  - ZZ_test[[k]] %*% (identified_ijm$alpha_0_hat  + identified_ijm$alpha_k_hat[[k]]))^2)/TT_test 
          }
          loss_eta_thres[i,j,m] <- tmp_sum_ijm
        }
      }
    }
    list_block_error[[l]] <- loss_eta_thres
    
    cat(l,"th fold done.\n")
  }
  
  # ------------------------------------------- #
  # Choose eta and threshold_0:
  mean_error <- array(NA,dim=c(grid_eta,grid_thres_0,grid_thres_k))
  for (i in 1:grid_eta){
    for (j in 1:grid_thres_0){
      for (m in 1:grid_thres_k){
        mean_error[i,j,m] <- sum(mapply(l=1:nfold,function(l)list_block_error[[l]][i,j,m]))/nfold
      }
    }
  }
  idx <- which(mean_error==min(mean_error),arr.ind=TRUE)
  eta_hat <- eta_grid[idx[1,1]]
  thres_0_hat <- thres_0_grid[idx[1,2]]
  thres_k_hat <- thres_k_grid[,idx[1,3]]
  
  # ---------------------------------------------------- #
  # Return output:
  output <- list(eta_grid = eta_grid, 
                 thres_0_grid = thres_0_grid,
                 thres_k_grid = thres_k_grid,
                 eta_hat=eta_hat,
                 thres_0_hat = thres_0_hat,
                 thres_k_hat = thres_k_hat,
                 mean_error=mean_error)
  return(output)
}