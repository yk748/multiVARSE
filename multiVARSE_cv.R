multiVARSE_cv <- function(X_t,K,d,p,TT,nfold=5){
  
  eta_grid <- seq(0.1,1,0.1)
  C_0 <- seq(0.01,1,length=10)
  N_min <- K*min(TT)
  thres_0_grid <- C_0*sqrt(log(d^2*p)/N_min)
  loss_eta_thres <- array(NA,dim=c(length(eta_grid),length(thres_0_grid)))
  
  # for (k in 1:K){
  #   X_t[[k]] <- t(scale(t(X_t[[k]]),center=TRUE,scale=TRUE)) 
  # }
  list_block_error <- list()
  for (l in 1:nfold){
    
    # Training and test data for fixed fold 
    debiased_train <- list(); 
    TT_train_vec <- vector("numeric",K)
    vec_YY_test <- ZZ_test <- list()
    for (k in 1:K){
      
      # Indexing:
      block_size <- round(TT[k]/nfold)
      if (l != nfold){
        idx_test <- (block_size*(l-1)+1):(block_size*l)
      }else{
        idx_test <- (block_size*(l-1)+1):TT[k]
      }
      idx_train <- c(1:TT[k])[-idx_test]
      
      # Estimating & debiasing with train data:
      X_t_train <- X_t[[k]][,idx_train]
      TT_train_vec[k] <- dim(X_t_train)[2]
      opt <- multiVARSE_optimizer(X_t_train,d,p,TT_train_vec[k])
      debiased_train[[k]] <- multiVARSE_debiasor(opt$beta_hat,
                                                 X_t_train,d,p,TT_train_vec[k])
    
      # Get YY and XX from test data:
      X_t_test <- X_t[[k]][,idx_test]
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
    }
    
    # Compute loss functions with fixed eta and thresholds
    for (i in 1:length(eta_grid)){
      for (j in 1:length(thres_0_grid)){
        identified_ij <- multiVARSE_identifier(debiased_train,
                                               eta_grid[i],thres_0_grid[j],
                                               K,d,p,TT_train_vec)
        tmp_sum_ij <- 0
        for (k in 1:K){
          TT_test <- length(vec_YY_test[[k]])
          tmp_sum_ij <- tmp_sum_ij + sum((vec_YY_test[[k]]  - ZZ_test[[k]] %*% (identified_ij$alpha_0_hat  + identified_ij$alpha_k_hat[[k]]))^2)/TT_test
        }
        loss_eta_thres[i,j] <- tmp_sum_ij
      }
    }
    list_block_error[[l]] <- loss_eta_thres
  }
  
  # Choose eta and threshold_0
  mean_error <- array(NA,dim=c(length(eta_grid),length(thres_0_grid)))
  for (i in 1:length(eta_grid)){
    for (j in 1:length(thres_0_grid)){
      mean_error[i,j] <- sum(mapply(l=1:nfold,function(l)list_block_error[[l]][i,j]))/nfold
    }
  }
  idx <- which(mean_error==min(mean_error),arr.ind=TRUE)
  eta_hat <- eta_grid[idx[1,1]]
  thres_hat <- thres_0_grid[idx[1,2]]
  
  output <- list(eta_hat=eta_hat,thres_hat=thres_hat,mean_error=mean_error)
  return(output)
}