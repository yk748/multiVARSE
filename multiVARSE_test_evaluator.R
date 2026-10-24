###################################################################################
# The following is the code for evaluating decisions in hypothesis tests 
# displayed in Section 5.3.
# The corresponding explanation is included in Section S.4.1 
# in the Supplementary Material
###################################################################################
multiVARSE_test_evaluator <- function(data,test,K,d,p,type,alpha){
  
  # ------------------------------------------- #
  # Compute FDR and power:
  ## Compute T0/T0c/TK/TKc:
  T0 <- which(data$A_0!=0,arr.ind=TRUE)
  T0c <- which(data$A_0==0,arr.ind=TRUE)
  TK <- c(); TKc <- c()
  for (k in 1:K){
    TK <- rbind(TK,which(data$A_k[[k]]!=0,arr.ind=TRUE))
    TKc <- rbind(TKc,which(data$A_k[[k]]==0,arr.ind=TRUE))
  }
  TK <- unique(TK); TKc <- unique(TKc)
  TKc <- TKc[!duplicated(rbind(TK, TKc))[-seq_len(nrow(TK))], ]
  
  # ------------------------------------------- #
  ## Compute T0K/T0cK/T0Kc/T0cKc:
  T0K <- intersect(split(T0, row(T0)), split(TK, row(TK)))
  T0K <- do.call(rbind, lapply(T0K, matrix, ncol = 2))
  
  T0Kc <- intersect(split(T0, row(T0)), split(TKc, row(TKc)))
  T0Kc <- do.call(rbind, lapply(T0Kc, matrix, ncol = 2))
  
  T0cK <- intersect(split(T0c, row(T0c)), split(TK, row(TK)))
  T0cK <- do.call(rbind, lapply(T0cK, matrix, ncol = 2))
  
  T0cKc <- intersect(split(T0c, row(T0c)), split(TKc, row(TKc)))
  T0cKc <- do.call(rbind, lapply(T0cKc, matrix, ncol = 2))
  
  if (type == "nullity"){
    true_alt <- rbind(T0K,T0cK,T0Kc)
    true_null <- T0cKc
  }else if (type == "homogeneity"){
    true_alt <- rbind(T0K,T0cK)
    true_null <- rbind(T0Kc,T0cKc)
  }else if (type == "significance"){
    true_alt <- rbind(T0K,T0Kc)
    true_null <- rbind(T0cK,T0cKc)
  }
  
  # Compute FDR and power:
  res_uncontrolled <- compute_fdr_power(test$reject_idx,
                                        true_alt, true_null, d, p)
  res_controlled   <- compute_fdr_power(test$reject_idx_controlled,
                                        true_alt, true_null, d, p)
  
  # ------------------------------------------- #
  # Return output:
  output <- list(
    type      = type,
    true_alt  = true_alt,
    true_null = true_null,
    hat_zeta  = test$hat_zeta,
    # Uncontrolled:
    FDR           = res_uncontrolled$FDR,
    power         = res_uncontrolled$power,
    # Controlled:
    FDR_controlled   = res_controlled$FDR,
    power_controlled = res_controlled$power
  )
  return(output)
}


#################################################################################
# Helper: compute FDR and power given a rejection index set
compute_fdr_power <- function(rej_idx, true_alt, true_null, d, p) {
  n_rej <- if (is.null(nrow(rej_idx)) || nrow(rej_idx) == 0) 0L else nrow(rej_idx)
  n_alt <- if (is.null(nrow(true_alt)) || nrow(true_alt) == 0) 1L else nrow(true_alt)
  
  true_mat <- dec_mat <- matrix(0, nrow = d, ncol = (d * p))
  true_mat[true_null] <- 1
  if (n_rej > 0) dec_mat[rej_idx] <- 1
  fdr <- sum(true_mat * dec_mat) / max(1, n_rej)
  
  true_mat <- dec_mat <- matrix(0, nrow = d, ncol = (d * p))
  true_mat[true_alt] <- 1
  if (n_rej > 0) dec_mat[rej_idx] <- 1
  pwr <- sum(true_mat * dec_mat) / n_alt
  
  list(FDR = fdr, power = pwr)
}
