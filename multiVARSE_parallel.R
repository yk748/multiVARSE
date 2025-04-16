
get_sim <- function(K,d,TT_vec,p,level,iter_max){
  
  if (level == "low"){
    s0 <- 0.05; sk <- 0.01
  }else if (level == "med") {
    s0 <- 0.03; sk <- 0.03
  }else if (level == "high"){
    s0 <- 0.01; sk <- 0.05 
  }
  
  # 10%
  
  library("doParallel")
  library("dqrng")
  library("doRNG")
  
  list_dummy <- list()
  cl <- makeCluster(detectCores())
  # cl <- makeCluster(2)
  registerDoParallel(cl)
  
  set.seed(123)
  multiVARSE_sim <- foreach(iter = 1:iter_max,
                      .errorhandling = 'remove',
                      .packages = c("mvtnorm","glmnet","Matrix","multivar")) %dopar% {
                        
                        source("multiVARSE_simulator.R")
                        source("multiVARSE_optimizer.R")
                        source("multiVARSE_debiasor.R")
                        source("multiVARSE_identifier.R")
                        source("multiVARSE_cv.R")
                        source("multiVARSE_tester.R")
                        source("multiVARSE_benchmark.R")  
                        
                        list_data <- multiVARSE_simulator(K,d,p,TT_vec,s0,sk,eps=1)
                        
                        list_opt <- list()
                        list_debiased <- list()
                        list_identified <- list()
                        for (k in 1:K){
                          # Individual estimation:
                          X_t <- list_data$x_tk[[k]]
                          list_opt[[k]] <- multiVARSE_optimizer(X_t,d,p,TT_vec[k])
                          
                          # Individual debiasing:
                          beta_hat <- list_opt[[k]]$beta_hat
                          list_debiased[[k]] <- multiVARSE_debiasor(beta_hat,X_t,d,p,TT_vec[k])
                        }
                        
                        # Cross validation for thresholds
                        cv_tuning <- multiVARSE_cv(list_data$x_tk,K,d,p,TT_vec,nfold=5)
                        
                        # Identifying common and unique components
                        eta_hat <- cv_tuning$eta_hat
                        thres_0_hat <- cv_tuning$thres_hat
                        list_identified <- multiVARSE_identifier(list_debiased,
                                                                 eta_hat,thres_0_hat,
                                                                 K,d,p,TT_vec)
                        
                        # Hypotheses tests on 
                        test_sig <- multiVARSE_tester(list_data,list_debiased,
                                                      type="significance",K,d,p,TT_vec)
                        test_hetero <- multiVARSE_tester(list_data,list_debiased,
                                                         type="heterogeneity",K,d,p,TT_vec)
                        
                        # Benchmark method: multivar
                        list_benchmark <- multiVARSE_benchmark(list_data,K,d)
                        
                        # Measure: Estimation
                        beta_k <- mapply(k=1:K,function(k){list_identified$alpha_0_hat + list_identified$alpha_k_hat[[k]]})
                        list_measure <- list()
                        
                        # RMSE
                        rmse <- list(alpha_0 =  norm(list_identified$alpha_0_hat - list_data$alpha_0,"2")/norm(list_data$alpha_0,"2"),
                        alpha_k = mapply(k=1:K,function(k){norm(list_identified$alpha_k_hat[[k]] 
                                                                - list_data$alpha_k[[k]],"2")/norm(list_data$alpha_k[[k]],"2")}),
                        beta_k = mapply(k=1:K,function(k)norm(beta_k[,k]  
                                                              - list_data$beta_k[[k]],"2")/norm(list_data$beta_k[[k]],"2")))
                        
                        rmse_bench <- list(alpha_0 =  norm(as.vector(t(list_benchmark$common)) 
                                             - list_data$alpha_0,"2")/norm(list_data$alpha_0,"2"),
                       alpha_k = mapply(k=1:K,function(k){norm(as.vector(t(list_benchmark$unique[[k]])) 
                                                               - list_data$alpha_k[[k]],"2")/norm(list_data$alpha_k[[k]],"2")}),
                       beta_k = mapply(k=1:K,function(k)norm(as.vector(t(list_benchmark$total[[k]]))
                                                             - list_data$beta_k[[k]],"2")/norm(list_data$beta_k[[k]],"2")))
                        
                        list_measure$rmse <- list(propose=rmse,benchmark=rmse_bench)
                        
                        # Sensitivity
                        sens <- list(alpha_0 = sum(list_identified$alpha_0_hat!=0 
                                                   & list_data$alpha_0!=0)/sum(list_data$alpha_0!=0),
                      alpha_k = mapply(k=1:K,function(k)sum(list_identified$alpha_k_hat[[k]]!=0 
                                                            & list_data$alpha_k[[k]]!=0)/sum(list_data$alpha_k[[k]]!=0)),
                      beta_k = mapply(k=1:K,function(k)sum(beta_k[,k]!=0 
                                                           & list_data$beta_k[[k]]!=0)/sum(list_data$beta_k[[k]]!=0)))
                        
                sens_bench <- list(alpha_0 = sum(as.vector(t(list_benchmark$common))!=0 
                                                 & list_data$alpha_0!=0)/sum(list_data$alpha_0!=0),
                   alpha_k = mapply(k=1:K,function(k)sum(as.vector(t(list_benchmark$unique[[k]]))!=0 
                                                         & list_data$alpha_k[[k]]!=0)/sum(list_data$alpha_k[[k]]!=0)),
                   beta_k = mapply(k=1:K,function(k)sum(as.vector(t(list_benchmark$total[[k]]))!=0 
                                                        & list_data$beta_k[[k]]!=0)/sum(list_data$beta_k[[k]]!=0)))
                        
                        list_measure$sens <- list(propose=sens,benchmark=sens_bench)
                        
                
                        # Specificity
                        spec <- list(alpha_0 = sum(list_identified$alpha_0_hat==0 
                                                   & list_data$alpha_0==0)/sum(list_data$alpha_0==0),
                        alpha_k = mapply(k=1:K,function(k)sum(list_identified$alpha_k_hat[[k]]==0 
                                                              & list_data$alpha_k[[k]]==0)/sum(list_data$alpha_k[[k]]==0)),
                        beta_k = mapply(k=1:K,function(k)sum(beta_k[,k]==0 
                                                             & list_data$beta_k[[k]]==0)/sum(list_data$beta_k[[k]]==0)))
                        
                spec_bench <- list(alpha_0 = sum(as.vector(t(list_benchmark$common))==0 
                                                 & list_data$alpha_0==0)/sum(list_data$alpha_0==0),
                                   alpha_k = mapply(k=1:K,function(k)sum(as.vector(t(list_benchmark$unique[[k]]))==0 
                                                                         & list_data$alpha_k[[k]]==0)/sum(list_data$alpha_k[[k]]==0)),
                                   beta_k = mapply(k=1:K,function(k)sum(as.vector(t(list_benchmark$total[[k]]))==0 
                                                                        & list_data$beta_k[[k]]==0)/sum(list_data$beta_k[[k]]==0)))
                
                        list_measure$spec <- list(propose=spec,benchmark=spec_bench)
                        
                        # Output
                        output <- list()
                        output$data <- list_data
                        output$opt <- list_opt
                        output$deb <- list_debiased
                        output$cv <- cv_tuning
                        output$idt <- list_identified
                        output$sig <- test_sig
                        output$hetero <- test_hetero
                        output$measure <- list_benchmark
                        output$bench <- list_benchmark
                        list_dummy[[iter]] <- output
                      }  
  stopCluster(cl)
  save(multiVARSE_sim,file=paste0("multiVARSE_d",d,"_T",TT_vec[1],"_level_",level,".RData"))
}

#####################################################################################
K <- 15; p <- 1
d_list <- c(10,30,50)
TT_list <- c(50,100,200)
hetero_list <- c("low","med","high")

for (i in 3:3){
  d <- d_list[i]
  
  for (j in 1:3){
    TT_vec <- rep(TT_list[j],K)
    
    for (k in 1:3){
      level <- hetero_list[k]
      
      cat("The current iteration is d:",d,"TT:",TT_list[j],"hetero level:",level,"\n")
      get_sim(K,d,TT_vec,p,level,iter_max=50)
      cat("Done.\n")

    } 
  }
}

#####################################################################################
K <- 15; p <- 1
d_list <- c(10,30,50)
TT_list <- c(50,100,200)
hetero_list <- c("low","med","high")

for (i in 3:3){
  d <- d_list[i]
  
  for (j in 1:3){
    TT_vec <- rep(TT_list[j],K)
    
    for (k in 1:3){
      level <- hetero_list[k]
      
      cat("The current iteration is d:",d,"TT:",TT_list[j],"hetero level:",level,"\n")
      load(paste0("multiVARSE_d",d,"_T",TT_vec[1],"_level_",level,".RData"))
    }
  }
}


# alpha_0, rmse
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)multiVARSE_sim[[iter]]$measure$rmse$propose$alpha_0))
hist(mapply(iter=1:50,function(iter)multiVARSE_sim[[iter]]$measure$rmse$benchmark$alpha_0))

# mean alpha_k, rmse
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$rmse$propose$alpha_k)))
     hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$rmse$benchmark$alpha_k)))

# mean beta_k, rmse
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$rmse$propose$beta_k)))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$rmse$benchmark$beta_k)))

# alpha_0, sens
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)multiVARSE_sim[[iter]]$measure$sens$propose$alpha_0))
hist(mapply(iter=1:50,function(iter)multiVARSE_sim[[iter]]$measure$sens$benchmark$alpha_0))

# mean alpha_k, sens
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$sens$propose$alpha_k)))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$sens$benchmark$alpha_k)))

# mean beta_k, sens
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$sens$propose$beta_k)))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$sens$benchmark$beta_k)))

# alpha_0, sens
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)multiVARSE_sim[[iter]]$measure$spec$propose$alpha_0))
hist(mapply(iter=1:50,function(iter)multiVARSE_sim[[iter]]$measure$spec$benchmark$alpha_0))

# mean alpha_k, sens
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$spec$propose$alpha_k)))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$spec$benchmark$alpha_k)))

# mean beta_k, sens
par(mfrow=c(2,1))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$spec$propose$beta_k)))
hist(mapply(iter=1:50,function(iter)mean(multiVARSE_sim[[iter]]$measure$spec$benchmark$beta_k)))




