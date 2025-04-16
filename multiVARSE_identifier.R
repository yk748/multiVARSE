multiVARSE_identifier <- function(debias,eta,thres_0,K,d,p,TT){
  
  N_min <- K*min(TT)
  thres_k <- vector("numeric",K)
  for (k in 1:K){
    thres_k[k] <- thres_0*sqrt(N_min/TT[k]) 
  }
  
  alpha_0_hat <- theta_tilde <- vector("numeric",d^2*p)
  for (j in 1:(d^2*p)){
    theta_j <- mapply(k=1:K,function(k)debias[[k]]$b_hat[j])
    
    grid <- seq(min(theta_j),max(theta_j),by=0.03)
    red_loss <- vector("numeric",length(grid))
    for (i in 1:length(grid)){
      red_loss[i] <- sum(mapply(k=1:K,function(k)red_ftn((theta_j[k]-grid[i]),eta)))
    }
    theta_tilde[j] <- grid[which.min(red_loss)]
    
    # Hard threshold
    alpha_0_hat[j] <- hard_thres(theta_tilde[j],thres_0)
    # alpha_0_hat[j] <- soft_thres(theta_tilde[j],thres_0)
  }
  
  alpha_k_hat <- delta_tilde <- list()
  for (k in 1:K){
    delta_tilde[[k]] <- debias[[k]]$b_hat - theta_tilde
    
    tmp_vec <- vector("numeric",d^2*p)
    for (j in 1:(d^2*p)){
      tmp_vec[j] <- hard_thres(delta_tilde[[k]][j],thres_k[k])
      # tmp_vec[j] <- soft_thres(delta_tilde[[k]][j],thres_k[k])
    }
    alpha_k_hat[[k]] <- tmp_vec
  }
  
  output <- list(theta_tilde = theta_tilde, delta_tilde = delta_tilde,
                 alpha_0_hat = alpha_0_hat, alpha_k_hat = alpha_k_hat,
                 thres_0 = thres_0, thres_k = thres_k)
  return(output)
}

#################################################################################
hard_thres <- function(val,thres){
  return(val*ifelse(abs(val)>=thres,1,0))
}

soft_thres <- function(val,thres){
  return(sign(val)*max(abs(val)-thres,0))
}

#################################################################################
objective_function <- function(c, theta, eta) {
  return(sum(sapply(theta - c, red_ftn, eta = eta)))
}

#################################################################################
# Huber:
huber_ftn <- function(x,eta){
  if (abs(x) <= eta){
    return(1/2*x^2)
  }else{
    return(eta*abs(x) - eta^2/2)
  }
}

# Tukey biweighted
tukey_ftn <- function(x,eta){
  if (abs(x) <= eta){
    return(eta^2/6*(1-(1-(x/eta)^2)^3))
  }else{
    return(eta^2/6)
  }
}
# redescending
red_ftn <- function(x,eta){
  return(min(x^2,eta^2))
}

# Lorentz
lor_ftn <- function(x,eta){
  return(log(1+x^2/eta))
}
# Welsch
welsch_ftn <- function(x,eta){
  return(eta^2/2*(1-exp(-x^2/eta^2)))
}
# ell1 + small quadratic
ql_ftn <- function(x,alpha){
  return(abs(x) + alpha*x^2)
}