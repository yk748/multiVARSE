multiVARSE_identifier <- function(debias,eta,thres_0,thres_k,thres_type,K,d,p,TT){
  
  # ------------------------------------------- #
  # Find alpha_0_tilde and alpha_0_hat:
  beta_mat <- sapply(debias, function(obj) obj$beta_tilde)  
  ind_contribute <- list()
  alpha_0_hat <- alpha_0_tilde <- vector("numeric",d^2*p)
  for (j in 1:(d^2*p)) {
    beta_tilde_j <- beta_mat[j,]
    
    grid <- seq(min(beta_tilde_j), max(beta_tilde_j), by = 0.01)
    red_loss <- colSums(outer(beta_tilde_j, grid, function(b,g) red_ftn(b-g, eta)))
    alpha_0_tilde[j] <- grid[which.min(red_loss)]
    ind_contribute[[j]] <- which(abs(beta_tilde_j - alpha_0_tilde[j]) <= eta)
    
    # thresholding
    if (thres_type == "hard") {
      alpha_0_hat[j] <- hard_thres(alpha_0_tilde[j], thres_0)
    } else {
      alpha_0_hat[j] <- soft_thres(alpha_0_tilde[j], thres_0)
    }
  }
  
  # ------------------------------------------- #
  # Find alpha_k_tilde and alpha_k_hat:
  alpha_k_hat <- alpha_k_tilde <- list()
  for (k in 1:K) {
    alpha_k_tilde[[k]] <- debias[[k]]$beta_tilde - alpha_0_tilde
    
    if (thres_type == "hard") {
      alpha_k_hat[[k]] <- hard_thres(alpha_k_tilde[[k]], thres_k[k])
    } else {
      alpha_k_hat[[k]] <- soft_thres(alpha_k_tilde[[k]], thres_k[k])
    }
  }
  
  # ---------------------------------------------------- #
  # Return output:
  output <- list(alpha_0_tilde = alpha_0_tilde, 
                 alpha_k_tilde = alpha_k_tilde,
                 ind_contribute = ind_contribute,
                 alpha_0_hat = alpha_0_hat, 
                 alpha_k_hat = alpha_k_hat,
                 thres_0 = thres_0, 
                 thres_k = thres_k)
  return(output)
}

#################################################################################
# Thresholds:
hard_thres <- function(val,thres){
  return(val*ifelse(abs(val)>=thres,1,0))
}

soft_thres <- function(val,thres){
  return(sign(val)*max(abs(val)-thres,0))
}

#################################################################################
# Get objective:
objective_function <- function(c, theta, eta) {
  return(sum(sapply(theta - c, red_ftn, eta = eta)))
}

red_ftn <- function(x, eta) {
  pmin(x^2, eta^2)
}
