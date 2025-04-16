multiVARSE_optimizer <- function(X_t,d,p,TT){
  
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
  vec_YY <- as.vector(YY)
  ZZ <- bdiag(replicate(d,XX,simplify=FALSE))
  
  cv.opt <- cv.glmnet(x=ZZ,y=vec_YY,intercept=FALSE,alpha=1)
  lambda <- cv.opt$lambda.min
  opt <- glmnet(x=ZZ,y=vec_YY,intercept=FALSE,alpha=1,lambda=lambda)
  beta_hat <- opt$beta
  Phi_hat <- array(unlist(lapply(1:p,function(x){t(array(beta_hat,dim=c(d,d,p))[,,x])})),dim=c(d,d,p))
  
  output <- list(Phi_hat = Phi_hat, beta_hat = beta_hat, lambda_cv = lambda)
  return(output)
}
