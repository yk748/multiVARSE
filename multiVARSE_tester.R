multiVARSE_tester <- function(data,debias,type,K,d,p,TT){
  
  # Set up:
  Cov_X <- list()
  N_k <- vector("numeric",K)
  V_hat <- list()
  for (k in 1:K){
    N_k[k] <- TT[k] - p
    Cov_X[[k]] <- data$x_tk[[k]] %*% t(data$x_tk[[k]])/TT[k]
    V_hat[[k]] <- debias[[k]]$Theta_hat %*% Cov_X[[k]] %*% debias[[k]]$Theta_hat
  }
  
  if (type == "significance"){
    D <- diag(1,K)
  }else if (type == "heterogeneity"){
    D <- matrix(0,(K-1),K)
    for (k in 1:(K-1)){
      D[k,k] <- 1
      D[k,(k+1)] <- -1
    } 
  }
  
  scale_mat <- list()
  if (type == "significance"){
    for (j in 1:d){
      scale_mat[[j]] <-  diag(1/sqrt(mapply(k=1:K,function(k)V_hat[[k]][j,j]/N_k[k])))
    }
  }else if (type == "heterogeneity"){
    for (j in 1:d){
      scale_mat[[j]] <-  sqrt(solve(D %*% diag(mapply(k=1:K,function(k)V_hat[[k]][j,j]/N_k[k])) %*% t(D)))
    }
  }
  
  S_j <- list()
  idx <- (1:(d^2*p) %% d)
  idx[which(idx==0)] <- d
  for (j in 1:(d^2*p)){
    S_j[[j]] <- scale_mat[[idx[j]]] %*%  D %*% mapply(k=1:K,function(k)debias[[k]]$b_hat[j])
  }
  
  if (type == "significance"){
    df <- K
    fill_tab <- array(NA,dim=c((d^2*p),K))
    for (k in 1:K){
      fill_tab[,k] <- data$beta_k[[k]]
    }
    
    alpha <- sum((rowSums(t(mapply(j=1:(d^2*p),function(j)D %*% fill_tab[j,]))) == 0) 
                 & (mapply(j=1:(d^2*p),function(j)norm(S_j[[j]],"2")^2) >= qchisq(0.95,df))) / sum(
                   rowSums(t(mapply(j=1:(d^2*p),function(j)D %*% fill_tab[j,]))) == 0)
    power <- sum((rowSums(t(mapply(j=1:(d^2*p),function(j)D %*% fill_tab[j,]))) != 0) 
                 & (mapply(j=1:(d^2*p),function(j)norm(S_j[[j]],"2")^2) >= qchisq(0.95,df))) / sum(
                   rowSums(t(mapply(j=1:(d^2*p),function(j)D %*% fill_tab[j,]))) != 0)
    
  }else if (type == "heterogeneity"){
    df <- K-1
    fill_tab <- array(NA,dim=c((d^2*p),K))
    for (k in 1:K){
      fill_tab[,k] <- data$alpha_k[[k]]
    }
    idx <- setdiff(1:(d^2*p),which(data$alpha_0!=0))
    
    alpha <- sum((rowSums(t(mapply(j=idx,function(j)D %*% fill_tab[j,]))) == 0) 
                 & (mapply(j=idx,function(j)norm(S_j[[j]],"2")^2) >= qchisq(0.95,df))) / sum(
                   rowSums(t(mapply(j=idx,function(j)D %*% fill_tab[j,]))) == 0)
    power <- sum((rowSums(t(mapply(j=idx,function(j)D %*% fill_tab[j,]))) != 0) 
                 & (mapply(j=idx,function(j)norm(S_j[[j]],"2")^2) >= qchisq(0.95,df))) / sum(
                   rowSums(t(mapply(j=idx,function(j)D %*% fill_tab[j,]))) != 0)
  }
  
  output <- list(S_j=S_j,scale_mat=scale_mat,alpha=alpha,power=power)
  return(output)
}