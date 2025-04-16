multiVARSE_simulator <- function(K,d,p,TT,s0,sk,eps=1){
  
  # ---------------------------------------------------- #
  # Assigning the entries/values of the common/unique components
  prop.comm <- s0
  prop.uniq <- sk
  num.comm <- round(d^2*p*prop.comm)
  num.uniq <- round(d^2*p*prop.uniq)
  
  # Currently, only the case when overlapping does not happen.
  num.tot <- num.comm + num.uniq
  
  idx.comm <- sample(seq(1:(d^2*p)),num.comm,replace=FALSE)
  entry.comm <- c(runif(num.comm,min=0.1,max=0.9),vector("numeric",length=(d^2*p-num.comm)))
  entry.comm <- entry.comm[order( c( idx.comm, seq(1:(d^2*p))[-c(idx.comm)] ) )]
  
  idx.uniq <- array(NA,dim=c(num.uniq,K))
  entry.uniq <- array(NA,dim=c((d^2*p),K))
  entry.total <- array(NA,dim=c((d^2*p),K))
  
  idx_cumm <- idx.comm
  for (k in 1:K){
    idx.uniq[,k] <- sample(seq(1:(d^2*p))[-idx_cumm],num.uniq,replace=FALSE)
    entry.uniq[,k] <- c(runif(num.uniq,min=0.1,max=0.9),vector("numeric",length=d^2*p-num.uniq))
    entry.uniq[,k] <- entry.uniq[order( c( idx.uniq[,k], seq(1:(d^2*p))[-c(idx.uniq[,k])] ) ),k]
    entry.total[,k] <- entry.comm + entry.uniq[,k]
    
    idx_cumm <- c(idx_cumm,idx.uniq[,k])
  }
  
  # ---------------------------------------------------- #
  # Scale down for all subjects to be stable
  test.mat <- array(entry.total,dim=c(d*p,d*p,K))
  test.eigen <- vector("numeric",K)
  for (k in 1:K){
    test.mat[,,k] <- companion_form_phi(array(entry.total[,k],dim=c(d,d,p)),d,p)
    test.eigen[k] <- max(abs(eigen(test.mat[,,k])$values))
  }
  test.eigen.max <- max(test.eigen)
  decay.rate <- 1
  if( test.eigen.max > 0.95 ){
    
    decay.rate <- 0.9
    while ( test.eigen.max > 0.95 ){
      eig.mat.comm <- companion_form_phi(decay.rate*array(entry.comm,dim=c(d,d,p)),d,p)
      eig.mat.uniq <- array(NA,dim=c(d*p,d*p,K))
      eig.mat.tot <- array(NA,dim=c(d*p,d*p,K))
      for (k in 1:K){
        eig.mat.uniq[,,k] <- companion_form_phi(decay.rate*array(entry.uniq[,k],dim=c(d,d,p)),d,p)
        eig.mat.tot[,,k] <- eig.mat.comm + eig.mat.uniq[,,k]
      }
      
      decay.rate <- decay.rate*0.9
      for (k in 1:K){
        test.eigen[k] <- max(abs(eigen(eig.mat.tot[,,k])$values))
      }
      test.eigen.max <- max(test.eigen)
    }
  }
  
  alpha_0 <- entry.comm*decay.rate
  A_0 <- array(unlist(lapply(1:p,
                             function(x){t(array(alpha_0,
                                                 dim=c(d,d,p))[,,x])})),
               dim=c(d,d,p))
  
  A_k <- list(); Phi_k <- list()
  alpha_k <- list(); beta_k <- list()
  for (k in 1:K){
    alpha_k[[k]] <- entry.uniq[,k]*decay.rate
    beta_k[[k]] <- alpha_0 + alpha_k[[k]]
    
    A_k[[k]] <- array(unlist(lapply(1:p,
                                    function(x){t(array(alpha_k[[k]],
                                                        dim=c(d,d,p))[,,x])})),
                      dim=c(d,d,p))
    Phi_k[[k]] <- A_0 + A_k[[k]]
  }
  
  # ---------------------------------------------------- #
  # Generate sample paths
  sigma_eps <- diag(eps,d)
  
  Burn <- 500
  X_tk <- list()
  for (k in 1:K){
    E_tk <- t(rmvnorm((TT[k]+Burn), mean = rep(0,d), sigma = sigma_eps, method="eigen"))
    
    X_tmp <- array(0,dim=c(d,(TT[k]+Burn)))
    for (t in (p+1):(TT[k]+Burn)){
      X_tmp[,t] <- rowSums(sapply(c(1:p),
                                  function(x){rowSums(Phi_k[[k]][,,x]%*%X_tmp[,(t-x)])})) + E_tk[,t]
    }
    X_tk[[k]] <- X_tmp[,(Burn+1):(TT[k]+Burn)]
  }
  
  # ---------------------------------------------------- #
  # Return output:
  output <- list(Phi_k = Phi_k, A_0 = A_0, A_k = A_k,
                 beta_k = beta_k, alpha_0 = alpha_0, alpha_k = alpha_k,
                 x_tk = X_tk)
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




# if (type == "F2022"){
#   
# }else if (type == "F2024"){
#   if (level == 1){
#     pi_I <- 1; pi_P <- 1/4
#   }else if (level == 2){
#     pi_I <- c(1,2/3,1/3); pi_P <- c(1/5,1/10,1/20);
#   }else if (level == 3){
#     pi_I <- c(1/3,2/3,1); pi_P <- c(1/5,1/10,1/20);
#   }
#   
#   idx <- 1:(d^2*p)
#   arr_list <- list(); subj_list <- list()
#   for (iter in 1:length(pi_P)){
#     idx_chosen <- sort(sample(idx,size=(d^2*p)*pi_P[iter],replace=FALSE))
#     arr_list[[iter]] <- mapply(x=1:length(idx_chosen),
#                                function(x){convert_index(idx_chosen[x],d,p)})
#     subj_list[[iter]] <- sort(sample(1:K,size=K*pi_I[iter],replace=FALSE))
#     
#     idx <- setdiff(idx,idx_chosen)
#   }
#   
#   
# }