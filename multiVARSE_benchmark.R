multiVARSE_benchmark <- function(data,K,d){
  
 list_multivar <- list()
 # list_multivar[["mat_com"]] <- matrix(data$A_0,d,d)
 for (k in 1:K){
   list_multivar[["data"]][[k]] <- t(data$x_tk[[k]])
   # list_multivar[["mat_ind_unique"]][[k]] <- matrix(data$A_k[[k]],d,d)
   # list_multivar[["mat_ind_final"]][[k]] <- list_multivar[["mat_com"]] + list_multivar[["mat_ind_unique"]][[k]]
 } 
 
 multivar_model <- multivar::constructModel(data=list_multivar$data)
 cv_multivar <- multivar::cv.multivar(multivar_model)
 
 return(cv_multivar$mats)
}







# # Problem size
# N_sum <- sum(TT)
# q <- d^2*p
# 
# # Matrix piling
# vec_YY_aug <- vector("numeric",N_sum*d*K)
# ZZ_aug <- Matrix(0,nrow=N_sum*d*K,ncol=q*(K+1),sparse=TRUE)
# for (k in 1:K){
#   X_t <- data$x_tk[[k]]
#   X_t_scaled <- X_t
#   
#   YY <- t(X_t_scaled[,TT[k]:(p+1)])
#   if ( p == 1 ){
#     XX <- t(X_t_scaled[,(TT[k]-1):p])
#   }else{
#     tmp_X <- t(X_t_scaled[,(TT[k]-1):p])
#     for (l in 2:p){
#       tmp_X <- cbind(tmp_X,t(X_t_scaled[,(TT[k]-l):(p-l)]))
#     }
#     XX <- tmp_X
#   }
#   vec_YY <- as.vector(YY)
#   ZZ <- bdiag(replicate(d,XX,simplify=FALSE))    
#   
#   y_idx <- ((k-1)*((TT[k]-p)*d)+1):(k*(TT[k]-p)*d)
#   x_idx <- ((k*q)+1):((k+1)*q)
#   vec_YY_aug[y_idx] <- vec_YY
#   ZZ_aug[y_idx,1:q] <- ZZ_aug[y_idx,x_idx] <- ZZ
# }
# 
# cv.opt <- cv.glmnet(x=ZZ_aug,y=vec_YY_aug,intercept=FALSE,alpha=1)
# lambda <- cv.opt$lambda.min
# opt <- glmnet(x=ZZ_aug,y=vec_YY_aug,intercept=FALSE,alpha=1,lambda=lambda)
# beta_hat <- opt$beta
# 
# 
# Phi_hat <- list()
# alpha0_hat <- array(unlist(lapply(1:p,function(x){t(array(beta_hat[1:q],dim=c(d,d,p))[,,x])})),dim=c(d,d,p))
# alphak_hat <- list()
# for (k in 1:K){
#   alphak_hat[[k]] <- beta_hat[(k*q):((k+1)*q)]
#   Phi_hat[[k]] <- array(unlist(lapply(1:p,function(x){t(array(alpha0_hat,dim=c(d,d,p))[,,x])})),dim=c(d,d,p)) + array(unlist(lapply(1:p,function(x){t(array(alphak_hat[[k]],dim=c(d,d,p))[,,x])})),dim=c(d,d,p))
# }
# 
# output <- list(Phi_hat = Phi_hat, beta_hat = beta_hat, lambda_cv = lambda)
# return(output)