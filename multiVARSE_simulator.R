multiVARSE_simulator <- function(K,d,p,TT,s0,sk,eps=1){
  
  # ---------------------------------------------------- #
  # Assigning the entries/values of the common/unique components:
  prop.comm <- s0
  prop.uniq <- sk
  
  # these should be list types.
  num.comm <- round(d^2*p*prop.comm)
  num.uniq <- num.tot <- list()
  for (k in 1:K){
    num.uniq[[k]] <- round(d^2*p*prop.uniq[k])
    num.tot[[k]] <-  num.comm + num.uniq[[k]]
  }
  
  # idx.uniq <- array(NA,dim=c(num.uniq,K))
  idx.uniq <- list()
  entry.uniq <- array(NA,dim=c((d^2*p),K))
  entry.total <- array(NA,dim=c((d^2*p),K))
  
  # ---------------------------------------------------- #
  # Picking up common entries + unique entries for the first subject:
  test.eigen <- 1
  while ( test.eigen > 0.95 ){
    # ---------------------------------------------------- #
    # Picking up common entries:
    idx.comm <- sample(seq(1:(d^2*p)),num.comm,replace=FALSE)
    entry.comm <- c(runif(num.comm,min=0.1,max=0.9),vector("numeric",length=(d^2*p-num.comm)))
    entry.comm <- entry.comm[order( c( idx.comm, seq(1:(d^2*p))[-c(idx.comm)] ) )]
    idx_cumm <- idx.comm
    
    # ---------------------------------------------------- #
    # Picking up unique entries for the first subject:
    idx.uniq[[1]] <- sample(seq(1:(d^2*p))[-idx_cumm],num.uniq[[1]],replace=FALSE)
    entry.uniq[,1] <- c(runif(num.uniq[[1]],min=0.1,max=0.9),
                        vector("numeric",length=d^2*p-num.uniq[[1]]))
    entry.uniq[,1] <- entry.uniq[order( c( idx.uniq[[1]], 
                                           seq(1:(d^2*p))[-c(idx.uniq[[1]])] ) ),1]
    entry.total[,1] <- entry.comm + entry.uniq[,1]
    idx_cumm <- c(idx_cumm,idx.uniq[[1]])
    
    eig.mat.comm <- companion_form_phi(array(entry.comm,dim=c(d,d,p)),d,p)
    
    eig.mat <- eig.mat.comm + companion_form_phi(array(entry.uniq[,1],dim=c(d,d,p)),d,p)
    test.eigen <- max(abs(eigen(eig.mat)$values))
  }
  cat("1 th subject done\n")
  
  # ---------------------------------------------------- #
  # Picking up unique entries for the rest of the subjects:
  for (k in 2:K){
    test.eigen <- 1
    while ( test.eigen > 0.95 ){
      # ---------------------------------------------------- #
      # Picking up unique entries:
      idx.uniq[[k]] <- sample(seq(1:(d^2*p))[-idx_cumm],num.uniq[[k]],replace=FALSE)
      entry.uniq[,k] <- c(runif(num.uniq[[k]],min=0.1,max=0.9),
                          vector("numeric",length=d^2*p-num.uniq[[k]]))
      entry.uniq[,k] <- entry.uniq[order( c( idx.uniq[[k]], 
                                             seq(1:(d^2*p))[-c(idx.uniq[[k]])] ) ),k]
      entry.total[,k] <- entry.comm + entry.uniq[,k]
      idx_cumm_tmp <- c(idx_cumm,idx.uniq[[k]])
      
      eig.mat <- eig.mat.comm + companion_form_phi(array(entry.uniq[,k],dim=c(d,d,p)),d,p)
      test.eigen <- max(abs(eigen(eig.mat)$values))
    }
    cat(k,"th subject done\n")
    idx_cumm <- idx_cumm_tmp
  }
  
  # ---------------------------------------------------- #
  # Define alpha0, alpha_k, and beta_k:
  alpha_k <- list(); beta_k <- list()
  
  alpha_0 <- entry.comm
  for (k in 1:K){
    alpha_k[[k]] <- entry.uniq[,k]
    beta_k[[k]] <- alpha_0 + alpha_k[[k]]
  }
  
  # ---------------------------------------------------- #
  # Define A_0,A_k, and Phi_k:
  A_k <- list(); Phi_k <- list()
  
  A_0 <- matrix(alpha_0,(d*p),d,byrow=FALSE)
  for (k in 1:K){
    A_k[[k]] <- matrix(alpha_k[[k]],(d*p),d,byrow=FALSE)
    Phi_tmp <- c()
    for (h in 1:p){
      Phi_tmp <- cbind(Phi_tmp,t(A_k[[k]][((h-1)*d+1):(h*d),] + A_0[((h-1)*d+1):(h*d),]))
    }
    Phi_k[[k]] <- array(Phi_tmp,c(d,d,p))
  }
  
  # ---------------------------------------------------- #
  # Generate sample paths:
  sigma_eps <- diag(eps,d)
  
  Burn <- 500
  X_t <- list()
  for (k in 1:K){
    E_tk <- t(rmvnorm((TT[k]+Burn), mean = rep(0,d), sigma = sigma_eps, method="eigen"))
    
    X_tmp <- array(0,dim=c(d,(TT[k]+Burn)))
    for (t in (p+1):(TT[k]+Burn)){
      X_tmp[,t] <- rowSums(sapply(c(1:p),
                                  function(x){rowSums(Phi_k[[k]][,,x]%*%X_tmp[,(t-x)])})) + E_tk[,t]
    }
    X_t[[k]] <- X_tmp[,(Burn+1):(TT[k]+Burn)]
  }
  
  # ---------------------------------------------------- #
  # Return output:
  output <- list(beta_k = beta_k, alpha_0 = alpha_0, alpha_k = alpha_k,
                 Phi_k = Phi_k, A_0 = A_0, A_k = A_k,
                 x_t = X_t)
  return(output)
}



#################################################################################
companion_form_phi = function(Psi,d,p){
  
  comp_Psi <- array(0,dim=c((p*d),(p*d)))
  if (p == 1){
    comp_Psi = Psi[,,1]
  }
  else{
    for (i in 1:p){
      for (j in 1:p){
        if(i == 1){
          comp_Psi[1:d,((j-1)*d+1):(j*d)] <- Psi[,,j]
        }
        else if ( (i-1) == j & i != 1 ){
          comp_Psi[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- diag(1,d)
        }
      }
    }
  }
  
  return(comp_Psi)
}

