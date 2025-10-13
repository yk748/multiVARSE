multiVARSE_test_evaluator <- function(data,test,K,d,p,type,alpha){
  
  # ------------------------------------------- #
  # Compute FDR and power:
  ## Compute S0/S0c/Sk/Skc:
  S0 <- which(data$A_0==0,arr.ind=TRUE)
  S0c <- which(data$A_0!=0,arr.ind=TRUE)
  Sk <- c(); Skc <- c()
  for (k in 1:K){
    Sk <- rbind(Sk,which(data$A_k[[k]]==0,arr.ind=TRUE))
    Skc <- rbind(Skc,which(data$A_k[[k]]!=0,arr.ind=TRUE))
  }
  Skc <- unique(Skc)
  Sk <- unique(Sk)
  Sk <- Sk[!duplicated(rbind(Skc, Sk))[-seq_len(nrow(Skc))], ]
  
  # ------------------------------------------- #
  ## Compute S0Sk/S0cSk/S0cSk/S0cSkc:
  common_S0Sk <- intersect(split(S0, row(S0)), split(Sk, row(Sk)))
  common_S0Sk <- do.call(rbind, lapply(common_S0Sk, matrix, ncol = 2))
  
  common_S0Skc <- intersect(split(S0, row(S0)), split(Skc, row(Skc)))
  common_S0Skc <- do.call(rbind, lapply(common_S0Skc, matrix, ncol = 2))
  
  common_S0cSk <- intersect(split(S0c, row(S0c)), split(Sk, row(Sk)))
  common_S0cSk <- do.call(rbind, lapply(common_S0cSk, matrix, ncol = 2))
  
  common_S0cSkc <- intersect(split(S0c, row(S0c)), split(Skc, row(Skc)))
  common_S0cSkc <- do.call(rbind, lapply(common_S0cSkc, matrix, ncol = 2))
  
  if (type == "nullity"){
    true_alt <- rbind(common_S0Skc,
                      common_S0cSk,
                      common_S0cSkc)
    true_null <- common_S0Sk
    
  }else if (type == "homogeneity"){
    true_alt <- rbind(common_S0Skc,
                      common_S0cSkc)
    true_null <- rbind(common_S0Sk,
                       common_S0cSk)
  }else if (type == "significance"){
    true_alt <- rbind(common_S0cSkc,
                      common_S0cSk)
    true_null <- rbind(common_S0Skc,
                       common_S0Sk)
  }
  
  # ------------------------------------------- #
  # FDR and power (without control)
  reject_idx <- test$reject_idx
  nonreject_idx <- test$nonreject_idx
  
  # FDR:
  if (is.null(dim(reject_idx))){
    FDR <- NA
  }else{
    true_mat <- dec_mat <- matrix(0,nrow=d,ncol=(d*p))
    true_mat[true_null] <- 1
    dec_mat[reject_idx] <- 1
    
    FDR <- sum(true_mat*dec_mat)/dim(reject_idx)[1]
  }
  
  # ------------------------------------------- #
  # Power:
  true_mat <- dec_mat <- matrix(0,nrow=d,ncol=(d*p))
  true_mat[true_alt] <- 1
  dec_mat[reject_idx] <- 1
  power <- sum(true_mat*dec_mat)/dim(true_alt)[1]
  
  # ------------------------------------------- #
  # FDR and power (with control)
  if (type == "nullity"){
    pval <- 1-pchisq(test$test_stat,K)
  }else if (type == "homogeneity"){
    pval <- 1-pchisq(test$test_stat,K-1)
  }else if (type == "significance"){
    pval <- 2*(1-pnorm(abs(test$test_stat)))
  }
  decision_ctrl <- matrix(fdr_procedure(pval,alpha=alpha,method="BH")$reject,d,(d*p))
  
  reject_ctrl <- which(decision_ctrl==TRUE,arr.ind=TRUE)
  nonreject_ctrl <- which(decision_ctrl==FALSE,arr.ind=TRUE)
  
  # FDR:
  if (is.null(dim(reject_ctrl))){
    FDR_ctrl <- NA
  }else{
    true_mat <- dec_mat <- matrix(0,nrow=d,ncol=(d*p))
    true_mat[true_null] <- 1
    dec_mat[reject_ctrl] <- 1
    
    FDR_ctrl <- sum(true_mat*dec_mat)/dim(reject_idx)[1]
  }
  
  # ------------------------------------------- #
  # Power:
  true_mat <- dec_mat <- matrix(0,nrow=d,ncol=(d*p))
  true_mat[true_alt] <- 1
  dec_mat[nonreject_ctrl] <- 1
  power_ctrl <- sum(true_mat*dec_mat)/dim(true_alt)[1]
  
  # ---------------------------------------------------- #
  # Return output:
  output <- list(type = type,
                 true_alt = true_alt,
                 true_null = true_null,
                 FDR = FDR,
                 power = power, 
                 FDR_ctrl = FDR_ctrl,
                 power_ctrl = power_ctrl )
  return(output)
  
}


#################################################################################
fdr_procedure <- function(pvals, alpha, method = c("BH","BY")) {
  
  # ---------------------------------------------------- #
  # setup:
  m <- length(pvals)
  ord <- order(pvals)
  pvals_ordered <- pvals[ord]
  
  # ---------------------------------------------------- #
  # Sequential procedure:
  if (method == "BH") {
    crit <- (1:m) * alpha / m
  } else if (method == "BY") {
    c_m <- sum(1 / (1:m))
    crit <- (1:m) * alpha / (m*c_m)
  }
  
  reject_idx <- which(pvals_ordered <= crit)
  reject <- rep(FALSE, m)
  if (length(reject_idx) > 0) {
    k <- max(reject_idx)
    reject[ord[1:k]] <- TRUE
  }

  # ---------------------------------------------------- #
  # Return output:
  out <- list()
  out$reject <- reject
  out$qvalues <- p.adjust(pvals, method = method)
  return(out)
}
