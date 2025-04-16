multiVARSE_debiasor <- function(beta_hat,X_t,d,p,TT){
  
  N <- TT - p
  beta_hat <- beta_hat
  
  # X_t_scaled <- t(scale(t(X_t),center=TRUE,scale=TRUE))
  X_t_scaled <- X_t
  
  YY <- t(X_t_scaled[,TT:(p+1)])
  if ( p == 1 ){
    XX <- t(X_t_scaled[,(TT-1):p])
  }else{
    tmp_X <- t(X_t_scaled[,(TT-1):p])
    for (l in 2:p){
      tmp_X <- cbind(tmp_X,t(X_t_scaled[,(TT-l):(p-l)]))
    }
    XX <- tmp_X
  }
  
  nodewise <- nodewise_regression(XX,d,p,TT)
  gamma_sqr_inv <- diag(1/nodewise$tau_sqr,d)
  Gamma_hat <- nodewise$Gamma_hat
  
  Theta_hat <- gamma_sqr_inv %*% Gamma_hat
  b_hat <- vector("numeric",d^2*p)
  for (i in 1:d){
    YY_i <- YY[,i]
    idx <- (d*p*(i-1)+1):(d*p*i)
    beta_i <- beta_hat[idx]
    
    b_hat[idx] <- beta_i + (Theta_hat/N) %*% t(XX) %*% (YY_i - XX %*% beta_i)
  }
  
  output <- list(b_hat = b_hat, Theta_hat = Theta_hat)
  return(output)
}

#################################################################################
nodewise_regression <- function(XX,d,p,TT){
  
  N <- TT - p
  Gamma_hat <- diag(1,d)
  tau_sqr <- vector("numeric",d)
  for (j in 1:d){
    XX_mj <- XX[,-j]
    XX_j <- XX[,j]
    cv.opt <- cv.glmnet(x=XX_mj,y=XX_j,intercept=FALSE,alpha=1)
    lambda <- cv.opt$lambda.min
    opt <- glmnet(x=XX_mj,y=XX_j,intercept=FALSE,alpha=1,lambda=lambda)
    
    Gamma_hat[-j,j] <- -as.vector(opt$beta)
    tau_sqr[j] <- sum((XX_j - XX_mj %*% opt$beta)^2)/N + 2*lambda*sum(abs(opt$beta))
  }
  
  output <- list(Gamma_hat = Gamma_hat, tau_sqr = tau_sqr)
  return(output)
}


