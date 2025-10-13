multiVARSE_optimizer <- function(X_t,d,p,TT){
  
  # ---------------------------------------------------- #
  # Stack X_t into regression:
  X_t_scaled <- t(scale(t(X_t)))
  
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
  vec_YY <- as.vector(YY)
  ZZ <- bdiag(replicate(d,XX,simplify=FALSE))
  
  # ---------------------------------------------------- #
  # Fit the model:
  cv.opt <- cv.glmnet(x=ZZ,y=vec_YY,intercept=FALSE,alpha=1)
  lambda <- cv.opt$lambda.min
  opt <- glmnet(x=ZZ,y=vec_YY,intercept=FALSE,alpha=1,lambda=lambda)
  beta_hat <- opt$beta
  rownames(beta_hat) <- NULL
  
  # ---------------------------------------------------- #
  # Return output:
  B_hat <- matrix(beta_hat,(d*p),d)
  Phi_hat <- array(NA,c(d,d,p))
  for (h in 1:p){
    Phi_hat[,,h] <- B_hat[((h-1)*d+1):(h*d),]
  }
  
  output <- list(Phi_hat = Phi_hat, beta_hat = beta_hat, lambda_cv = lambda)
  return(output)
}
