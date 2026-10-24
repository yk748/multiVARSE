###################################################################################
# The following is the code for generating synthetic data used in 
# the results in Section S.4.3 in the Supplementary Material.
###################################################################################
multiVARSE_simulator_overlap <- function(K,d,p,TT,case,eta=0.5,eps=0.5){
  
  # ---------------------------------------------------- #
  # set up:
  if (case == "case1"){
    pi_subj <- c(1,1/5)
    pi_prop <- c(1/20,1/20)
    
  }else if (case == "case2"){
    pi_subj <- c(1,2/5)
    pi_prop <- c(1/20,1/20)
    
  }else if (case == "case3"){
    pi_subj <- c(1,3/5)
    pi_prop <- c(1/20,1/20)
    
  }
  
  G <- length(pi_prop)
  n_entries <- d^2*p
  n_per_trl <- floor(n_entries*pi_prop)  # number of entries per group
  
  # Check feasibility: total entries across groups must not exceed n_entries
  if (sum(n_per_trl) > n_entries) {
    stop("Total proportion of entries across groups exceeds 1. Reduce pi_prop.")
  }
  
  # ---------------------------------------------------- #
  # Rejection sampler: resample until all subjects are stationary
  test.eigen <- 1
  while (test.eigen > 0.95) {
    
    # ---------------------------------------------------- #
    # Sample non-overlapping entry indices across groups:
    # Permute all entry indices once, then assign sequentially to groups
    permuted_entries <- sample(seq(1:n_entries))
    list_entry   <- list()
    list_subject <- list()
    idx <- 1
    for (g in 1:G) {
      list_entry[[g]]   <- permuted_entries[idx:(idx + n_per_trl[g]-1)]
      list_subject[[g]] <- sample(seq(1:K), floor(K*pi_subj[g]), replace = FALSE)
      idx <- idx + n_per_trl[g]
    }
    
    # ---------------------------------------------------- #
    # Fill entry.total using group membership:
    entry.total <- array(0, dim = c(n_entries, K))
    which_list  <- lapply(1:K, function(x) {which(sapply(list_subject, function(y) x %in% y))})
    
    test.eigen_vec <- rep(0, K)
    for (k in 1:K) {
      run <- which_list[[k]]
      for (g in run) {
        entry.total[list_entry[[g]], k] <- runif(n_per_trl[g],min = 0.1,max = 0.9)
      }
      test.mat          <- companion_form_phi(array(entry.total[,k], dim=c(d,d,p)),d,p)
      test.eigen_vec[k] <- max(abs(eigen(test.mat)$values))
    }
    test.eigen <- max(test.eigen_vec)
    cat("Spectral radius:", test.eigen, "\n")
  }

  # ---------------------------------------------------- #
  # Identifying common/unique components
  entry.comm <- rep(0,(d^2*p))
  entry.uniq <- array(0,dim=c((d^2*p),K))
  for (j in 1:(d^2*p)){
    beta_j <- entry.total[j,]
    nonzeros <- which(beta_j!=0)
    if (length(nonzeros) == 0){
      entry.comm[j] <- 0
      for (k in 1:K){
        entry.uniq[j,k] <- 0
      }
    }else{
      grid <- seq(min(beta_j[nonzeros]),max(beta_j[nonzeros]),by=0.01)
      red_loss <- vector("numeric",length(grid))
      for (i in 1:length(grid)){
        red_loss[i] <- sum(mapply(k=1:K,function(k)red_ftn((beta_j[k]-grid[i]),eta)))
      }
      entry.comm[j] <- grid[which.min(red_loss)]
      for (k in 1:K){
        entry.uniq[j,k] <- entry.total[j,k] - entry.comm[j]
      }
    }
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
    # E_tk <- matrix(rlaplace_clipped((TT[k]+Burn)*d, mu=0, b=1, lower=-5, upper=5),nrow=d)
    
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

#################################################################################
rlaplace_clipped <- function(n, mu, b, lower, upper) {
  
  u <- runif(n,-0.5,0.5)
  x <- mu - b*sign(u)*log(1-2*abs(u))
  
  # Step 2: clip to [lower, upper]
  x_clipped <- pmin(pmax(x, lower), upper)
  
  return(x_clipped)
}
