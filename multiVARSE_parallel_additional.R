###################################################################################
# The following is the loop for the simulation results displayed in Section S.4.3
# in the Supplementary Material.
###################################################################################
get_sim <- function(p,K=10,d,TT_mean,case,iter_max=50){
  
  set.seed(123)
  if (TT_mean == 100){
    TT_vec <- sample(seq(90,110),K,replace=TRUE)
  }else if (TT_mean == 200){
    TT_vec <- sample(seq(190,210),K,replace=TRUE)
  }

  library("doParallel")
  library("dqrng")
  library("doRNG")
  
  list_dummy <- list()
  available_cores <- detectCores()
  num_cores <- min(available_cores - 1, 50)
  num_cores <- max(num_cores, 1)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  sim_multivar <- foreach(iter = 1:iter_max,
                          .errorhandling = 'remove',
                          .packages = c("mvtnorm","glmnet","Matrix","multivar"),
                          .options.RNG = 123) %dopar% {
                            
                            source("multiVARSE_simulator_overlap.R")
                            source("multiVARSE_optimizer.R")
                            source("multiVARSE_debiasor.R")
                            source("multiVARSE_cv.R")
                            source("multiVARSE_identifier.R")
                            source("multiVARSE_benchmark.R")
                            source("multiVARSE_tester.R")
                            source("multiVARSE_test_evaluator.R")
                            
                            # --------------------------------------------------------------- #
                            # Data generation:
                            list_data <- multiVARSE_simulator_overlap(K,d,p,TT_vec,case=case,eta=0.5,eps=0.5)
                            
                            # --------------------------------------------------------------- #
                            # Benchmark methods: multiVAR & multiVAR (A)
                            data_multivar <- list()
                            for (k in 1:K){
                              data_multivar[["data"]][[k]] <- t(list_data$x_t[[k]])
                            }
                            
                            start_ls <- proc.time()
                            list_benchmark_ls <- multiVARSE_benchmark(data_multivar,K,d,type="non-adaptive")
                            time_ls <- proc.time() - start_ls
                            
                            start_als <- proc.time()
                            list_benchmark_als <- multiVARSE_benchmark(data_multivar,K,d,type="adaptive")
                            time_als <- proc.time() - start_als
                            
                            # --------------------------------------------------------------- #
                            # Evaluation - Estimation (Stratified Lasso):
                            rmse_alpha0_ls <- norm(as.vector(t(list_benchmark_ls$common))
                                                   - list_data$alpha_0,"2")/norm(list_data$alpha_0,"2")
                            rmse_mean_alpha_k_ls <- mean(mapply(k=1:K,function(k)
                              norm(as.vector(t(list_benchmark_ls$unique[[k]]))
                                   - list_data$alpha_k[[k]],"2")/norm(list_data$alpha_k[[k]],"2")))
                            rmse_mean_beta_k_ls <- mean(mapply(k=1:K,function(k)
                              norm(as.vector(t(list_benchmark_ls$total[[k]]))
                                   - list_data$beta_k[[k]],"2")/norm(list_data$beta_k[[k]],"2")))
                            
                            sens_mean_alpha_0_ls <- sum(as.vector(t(list_benchmark_ls$common))!=0
                                                        & list_data$alpha_0!=0)/sum(list_data$alpha_0!=0)
                            sens_mean_alpha_k_ls <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_ls$unique[[k]]))!=0
                                  & list_data$alpha_k[[k]]!=0)/sum(list_data$alpha_k[[k]]!=0)))
                            sens_mean_beta_k_ls <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_ls$total[[k]]))!=0
                                  & list_data$beta_k[[k]]!=0)/sum(list_data$beta_k[[k]]!=0)))
                            
                            spec_mean_alpha_0_ls <- sum(as.vector(t(list_benchmark_ls$common))==0
                                                        & list_data$alpha_0==0)/sum(list_data$alpha_0==0)
                            spec_mean_alpha_k_ls <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_ls$unique[[k]]))==0
                                  & list_data$alpha_k[[k]]==0)/sum(list_data$alpha_k[[k]]==0)))
                            spec_mean_beta_k_ls <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_ls$total[[k]]))==0
                                  & list_data$beta_k[[k]]==0)/sum(list_data$beta_k[[k]]==0)))
                            
                            
                            # Evaluation - Estimation (Adaptive Lasso):
                            rmse_alpha0_als <- norm(as.vector(t(list_benchmark_als$common))
                                                    - list_data$alpha_0,"2")/norm(list_data$alpha_0,"2")
                            rmse_mean_alpha_k_als <- mean(mapply(k=1:K,function(k)
                              norm(as.vector(t(list_benchmark_als$unique[[k]]))
                                   - list_data$alpha_k[[k]],"2")/norm(list_data$alpha_k[[k]],"2")))
                            rmse_mean_beta_k_als <- mean(mapply(k=1:K,function(k)
                              norm(as.vector(t(list_benchmark_als$total[[k]]))
                                   - list_data$beta_k[[k]],"2")/norm(list_data$beta_k[[k]],"2")))
                            
                            sens_mean_alpha_0_als <- sum(as.vector(t(list_benchmark_als$common))!=0
                                                         & list_data$alpha_0!=0)/sum(list_data$alpha_0!=0)
                            sens_mean_alpha_k_als <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_als$unique[[k]]))!=0
                                  & list_data$alpha_k[[k]]!=0)/sum(list_data$alpha_k[[k]]!=0)))
                            sens_mean_beta_k_als <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_als$total[[k]]))!=0
                                  & list_data$beta_k[[k]]!=0)/sum(list_data$beta_k[[k]]!=0)))
                            
                            spec_mean_alpha_0_als <- sum(as.vector(t(list_benchmark_als$common))==0
                                                         & list_data$alpha_0==0)/sum(list_data$alpha_0==0)
                            spec_mean_alpha_k_als <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_als$unique[[k]]))==0
                                  & list_data$alpha_k[[k]]==0)/sum(list_data$alpha_k[[k]]==0)))
                            spec_mean_beta_k_als <- mean(mapply(k=1:K,function(k)
                              sum(as.vector(t(list_benchmark_als$total[[k]]))==0
                                  & list_data$beta_k[[k]]==0)/sum(list_data$beta_k[[k]]==0)))
                            
                            # --------------------------------------------------------------- #
                            # Estimation:
                            
                            start_method <- proc.time()
                            
                            list_opt <- list()
                            list_debiased <- list()
                            list_identified <- list()
                            for (k in 1:K){
                              # Individual estimation:
                              X_t <- list_data$x_t[[k]]
                              list_opt[[k]] <- multiVARSE_optimizer(X_t,d,p,TT_vec[k])
                              
                              # Individual debiasing:
                              beta_hat <- list_opt[[k]]$beta_hat
                              list_debiased[[k]] <- multiVARSE_debiasor(beta_hat,X_t,d,p,TT_vec[k])
                              cat(k,"th subject done\n")
                            }
                            
                            # Cross validation for thresholds
                            cv_tuning <- multiVARSE_cv(list_data$x_t,list_debiased,list_opt,
                                                        agg_type = "redescending",
                                                        K,d,p,TT_vec,thres_type="hard",nfold=5)
                            
                            # Identifying common and unique components
                            eta_hat <- cv_tuning$eta_hat
                            thres_0_hat <- cv_tuning$thres_0_hat
                            thres_k_hat <- cv_tuning$thres_k_hat
                            
                            list_identified <- multiVARSE_identifier(list_debiased,
                                                                      agg_type = "redescending",
                                                                      eta_hat,thres_0_hat,thres_k_hat,
                                                                      thres_type="hard",K,d,p,TT_vec)
                            time_method <- proc.time() - start_method
                            
                            # --------------------------------------------------------------- #
                            # Evaluation - Estimation:
                            rmse_alpha0 <- norm(list_identified$alpha_0_hat
                                                - list_data$alpha_0,"2")/norm(list_data$alpha_0,"2")
                            rmse_mean_alpha_k <- mean(mapply(k=1:K,function(k)
                              norm(list_identified$alpha_k_hat[[k]]
                                   - list_data$alpha_k[[k]],"2")/norm(list_data$alpha_k[[k]],"2")))
                            rmse_mean_beta_k <- mean(mapply(k=1:K,function(k)
                              norm(list_identified$alpha_0_hat+list_identified$alpha_k_hat[[k]]
                                   - list_data$beta_k[[k]],"2")/norm(list_data$beta_k[[k]],"2")))
                            
                            sens_mean_alpha_0 <- sum(list_identified$alpha_0_hat!=0
                                                     & list_data$alpha_0!=0)/sum(list_data$alpha_0!=0)
                            spec_mean_alpha_0 <- sum(list_identified$alpha_0_hat==0
                                                     & list_data$alpha_0==0)/sum(list_data$alpha_0==0)
                            
                            sens_mean_alpha_k <- mean(mapply(k=1:K,function(k)sum(
                              (list_identified$alpha_k_hat[[k]])!=0
                              & list_data$alpha_k[[k]]!=0)/sum(list_data$alpha_k[[k]]!=0)))
                            spec_mean_alpha_k <- mean(mapply(k=1:K,function(k)sum(
                              (list_identified$alpha_k_hat[[k]])==0
                              & list_data$alpha_k[[k]]==0)/sum(list_data$alpha_k[[k]]==0)))
                            
                            sens_mean_beta_k <- mean(mapply(k=1:K,function(k)sum(
                              (list_identified$alpha_0_hat+list_identified$alpha_k_hat[[k]])!=0
                              & list_data$beta_k[[k]]!=0)/sum(list_data$beta_k[[k]]!=0)))
                            spec_mean_beta_k <- mean(mapply(k=1:K,function(k)sum(
                              (list_identified$alpha_0_hat+list_identified$alpha_k_hat[[k]])==0
                              & list_data$beta_k[[k]]==0)/sum(list_data$beta_k[[k]]==0)))
                            
                            # --------------------------------------------------------------- #
                            # Hypotheses tests:
                            test_sig <- multiVARSE_tester(list_data$x_t,list_debiased,
                                                          list_opt,list_identified,
                                                          type="significance",alpha=0.05,K,d,p,TT_vec)
                            test_homo <- multiVARSE_tester(list_data$x_t,list_debiased,
                                                           list_opt,list_identified,
                                                           type="homogeneity",alpha=0.05,K,d,p,TT_vec)
                            test_null <- multiVARSE_tester(list_data$x_t,list_debiased,
                                                           list_opt,list_identified,
                                                           type="nullity",alpha=0.05,K,d,p,TT_vec)
                            
                            # --------------------------------------------------------------- #
                            # Evaluation - Hypothesis tests:
                            eval_sig <- multiVARSE_test_evaluator(list_data,test_sig,K,d,p,
                                                                  type="significance",alpha=0.05)
                            eval_homo <- multiVARSE_test_evaluator(list_data,test_homo,K,d,p,
                                                                   type="homogeneity",alpha=0.05)
                            eval_null <- multiVARSE_test_evaluator(list_data,test_null,K,d,p,
                                                                   type="nullity",alpha=0.05)
                            
                            # --------------------------------------------------------------- #
                            # Output
                            output <- list()
                            output$data <- list_data
                            output$time <- list(ls=time_ls,als=time_als,method=time_method)
                            
                            # multiVAR:
                            output$est_bm_ls <- list_benchmark_ls
                            output$eval_est_ls <- c(rmse_alpha0_ls,
                                                    rmse_mean_alpha_k_ls,
                                                    rmse_mean_beta_k_ls,
                                                    sens_mean_alpha_0_ls,
                                                    sens_mean_alpha_k_ls,
                                                    sens_mean_beta_k_ls,
                                                    spec_mean_alpha_0_ls,
                                                    spec_mean_alpha_k_ls,
                                                    spec_mean_beta_k_ls)
                            
                            # multiVAR (A):
                            output$est_bm_als <- list_benchmark_als
                            output$eval_est_als <- c(rmse_alpha0_als,
                                                     rmse_mean_alpha_k_als,
                                                     rmse_mean_beta_k_als,
                                                     sens_mean_alpha_0_als,
                                                     sens_mean_alpha_k_als,
                                                     sens_mean_beta_k_als,
                                                     spec_mean_alpha_0_als,
                                                     spec_mean_alpha_k_als,
                                                     spec_mean_beta_k_als)
                            
                            output$opt <- list_opt
                            output$deb <- list_debiased
                            output$cv <- cv_tuning
                            output$idt <- list_identified
                            output$eval_est <- c(rmse_alpha0,
                                                 rmse_mean_alpha_k,
                                                 rmse_mean_beta_k,
                                                 sens_mean_alpha_0,
                                                 sens_mean_alpha_k,
                                                 sens_mean_beta_k,
                                                 spec_mean_alpha_0,
                                                 spec_mean_alpha_k,
                                                 spec_mean_beta_k)
                            
                            output$test_null <- test_null
                            output$test_homo <- test_homo
                            output$test_sig <- test_sig
                            output$eval_test <- c(eval_null$FDR_controlled,
                                                  eval_homo$FDR_controlled,
                                                  eval_sig$FDR_controlled,
                                                  eval_null$power_controlled,
                                                  eval_homo$power_controlled,
                                                  eval_sig$power_controlled)
                            list_dummy[[iter]] <- output
                          }  
  stopCluster(cl)
  print(length(sim_multivar))
  if (length(sim_multivar) == 0) {
    warning("No successful iterations were returned.")
  }
  save(sim_multivar,
       file=paste0("additional_multivarse_d",d,"_T",TT_mean,"_subject",K,"_",case,".RData"))
}

#####################################################################################
# Run simulations:
p <- 1

d_list <- c(10,20)
TT_list <- c(100,200)
case_list <- c("case1","case2","case3")

for (i in 1:2){
  d <- d_list[i]
  
  for (j in 1:2){
    TT_mean <- TT_list[j]
      
    for (k in 1:3){
      case <- case_list[k]
      
      cat("The current iteration is d:",d,",TT:",TT_mean,",subject: 10, case:",case,"\n")
      get_sim(p,K=10,d,TT_mean,case,iter_max=50)
      cat("Done.\n")
    } 
  }
}


###################################################################################
# The following is for generating Figure S.3 in Section S.4.3 
# in the Supplementary Material.
###################################################################################
library(ggplot2)
library(latex2exp)
library(patchwork)
library(Matrix)
library(RColorBrewer)

# --------------------------------------------------------------- #
# Estimation:
p <- 1

d_list <- c(10,20)
TT_list <- c(100,200)
case_list <- c("case1","case2","case3")
iter_max <- 50

df_est_prop <- data.frame()
df_est_ls <- data.frame()
df_est_als <- data.frame()

for (i in 1:2){
  d <- d_list[i]
  
  for (j in 1:2){
    TT_mean <- TT_list[j]
    
    for (k in 1:3){
      case <- case_list[k]

      cat("The current iteration is d:",d,",TT:",TT_mean,",subject: 10, case:",case,"\n")
      load(paste0("additional_multivarse_d",d,"_T",TT_mean,"_subject",10,"_",case,".RData"))
  
      # --------------------------------------------- #
      df_est_prop_tmp <- data.frame(t(mapply(iter=1:iter_max,function(iter)
        sim_multivar[[iter]]$eval_est)))
      colnames(df_est_prop_tmp) <- c("rmse_comm","rmse_uniq","rmse_tot",
                                     "sens_comm","sens_uniq","sens_tot",
                                     "spec_comm","spec_uniq","spec_tot")
      df_est_prop <- rbind(df_est_prop,
                           data.frame(df_est_prop_tmp,
                                      d = rep(d,iter_max),
                                      TT = rep(TT_mean,iter_max),
                                      case = rep(case,iter_max),
                                      method = rep("CE",iter_max)))
  
      # --------------------------------------------- #
      df_est_ls_tmp <- data.frame(t(mapply(iter=1:iter_max,function(iter)
        sim_multivar[[iter]]$eval_est_ls)))
      names(df_est_ls_tmp) <- c("rmse_comm","rmse_uniq","rmse_tot",
                                "sens_comm","sens_uniq","sens_tot",
                                "spec_comm","spec_uniq","spec_tot")
      df_est_ls <- rbind(df_est_ls,
                         data.frame(df_est_ls_tmp,
                                    d = rep(d,iter_max),
                                    TT = rep(TT_mean,iter_max),
                                    case = rep(case,iter_max),
                                    method = rep("ls",iter_max)))
  
      # --------------------------------------------- #
      df_est_als_tmp <- data.frame(t(mapply(iter=1:iter_max,function(iter)
        sim_multivar[[iter]]$eval_est_als)))
      names(df_est_als_tmp) <- c("rmse_comm","rmse_uniq","rmse_tot",
                                 "sens_comm","sens_uniq","sens_tot",
                                 "spec_comm","spec_uniq","spec_tot")
      df_est_als <- rbind(df_est_als,
                          data.frame(df_est_als_tmp,
                                     d = rep(d,iter_max),
                                     TT = rep(TT_mean,iter_max),
                                     case = rep(case,iter_max),
                                     method = rep("als",iter_max)))
    }
  }
}

# --------------------------------------------------------------- #
df_est <- rbind(df_est_prop,df_est_ls,df_est_als)
df_est$d <- factor(df_est$d,levels=unique(df_est$d),
                    labels=c("d=10","d=20"))
df_est$TT <- factor(df_est$TT,levels=unique(df_est$TT),
                    labels=c("T=100","T=200"))
df_est$case <- factor(df_est$case,levels=unique(df_est$case),
                      labels=c("Case1","Case2","Case3"))
df_est$method <- factor(df_est$method,levels=unique(df_est$method),
                        labels=c("Multi-VARSE","Multi-VAR","Multi-VAR (A)"))


df_rmse_comm <- df_est[,c(1,10,11,12,13)]
df_rmse_comm$comp <- factor(rep("Common",1800))
colnames(df_rmse_comm)[1] <- "rmse"
df_rmse_tot <- df_est[,c(3,10,11,12,13)]
df_rmse_tot$comp <- factor(rep("Individual",1800))
colnames(df_rmse_tot)[1] <- "rmse"
df_rmse <- rbind(df_rmse_comm,df_rmse_tot)

df_sens_comm <- df_est[,c(4,10,11,12,13)]
df_sens_comm$comp <- factor(rep("Common",1800))
colnames(df_sens_comm)[1] <- "sens"
df_sens_tot <- df_est[,c(6,10,11,12,13)]
df_sens_tot$comp <- factor(rep("Individual",1800))
colnames(df_sens_tot)[1] <- "sens"
df_sens <- rbind(df_sens_comm,df_sens_tot)

df_spec_comm <- df_est[,c(7,10,11,12,13)]
df_spec_comm$comp <- factor(rep("Common",1800))
colnames(df_spec_comm)[1] <- "spec"
df_spec_tot <- df_est[,c(9,10,11,12,13)]
df_spec_tot$comp <- factor(rep("Individual",1800))
colnames(df_spec_tot)[1] <- "spec"
df_spec <- rbind(df_spec_comm,df_spec_tot)


# RMSE:
p1 <- ggplot(df_rmse, aes(x = case, y = rmse, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ comp, scales = "free") +
  theme_minimal() +
  labs(y = TeX("RMSE"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Sensitivity:
p2 <- ggplot(df_sens, aes(x = case, y = sens, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ comp, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Sensitivity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom")

# Specificity:
p3 <- ggplot(df_spec, aes(x = case, y = spec, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ comp, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Specificity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom")


# --------------------------------------------------------------- #
# Hypothesis tests:
p <- 1

d_list <- c(10,20)
TT_list <- c(100,200)
case_list <- c("case1","case2","case3")
iter_max <- 50

df_test <- data.frame()

for (i in 1:2){
  d <- d_list[i]
  
  for (j in 1:2){
    TT_mean <- TT_list[j]
    
    for (k in 1:3){
      case <- case_list[k]
      
      cat("The current iteration is d:",d,",TT:",TT_mean,",subject: 10, case:",case,"\n")
      load(paste0("additional_multivarse_d",d,"_T",TT_mean,"_subject",10,"_",case,".RData"))
      
      # --------------------------------------------- #
      df_test_tmp <- data.frame(t(mapply(iter=1:iter_max,function(iter)
        sim_multivar[[iter]]$eval_test)))
      names(df_test_tmp) <- c("FDR_null","FDR_homo","FDR_sig",
                              "Power_null","Power_homo","Power_sig")
      df_test <- rbind(df_test,
                       data.frame(df_test_tmp,
                                  d = rep(d,iter_max),
                                  TT = rep(TT_mean,iter_max),
                                  case = rep(case,iter_max)))
      
    }
  }
}

# --------------------------------------------------------------- #
df_test$d <- factor(df_test$d,levels=unique(df_test$d),
                    labels=c("d=10","d=20"))
df_test$TT <- factor(df_test$TT,levels=unique(df_test$TT),
                     labels=c("T=100","T=200"))
df_est$case <- factor(df_est$case,levels=unique(df_est$case),
                      labels=c("Case1","Case2","Case3"))


df_fdr_sig <- df_test[,c(3,7,8,9)]
df_fdr_sig$comp <- factor(rep("FDR",600))
colnames(df_fdr_sig)[1] <- "perc"
df_power_sig <- df_test[,c(6,7,8,9)]
df_power_sig$comp <- factor(rep("Power",600))
colnames(df_power_sig)[1] <- "perc"
df_sig <- rbind(df_fdr_sig,df_power_sig)

df_fdr_homo <- df_test[,c(2,7,8,9)]
df_fdr_homo$comp <- factor(rep("FDR",600))
colnames(df_fdr_homo)[1] <- "perc"
df_power_homo <- df_test[,c(5,7,8,9)]
df_power_homo$comp <- factor(rep("Power",600))
colnames(df_power_homo)[1] <- "perc"
df_homo <- rbind(df_fdr_homo,df_power_homo)

df_fdr_null <- df_test[,c(1,7,8,9)]
df_fdr_null$comp <- factor(rep("FDR",600))
colnames(df_fdr_null)[1] <- "perc"
df_power_null <- df_test[,c(4,7,8,9)]
df_power_null$comp <- factor(rep("Power",600))
colnames(df_power_null)[1] <- "perc"
df_null <- rbind(df_fdr_null,df_power_null)


# Significance:
p4 <- ggplot(df_sig, aes(x = case, y = perc, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ comp, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Test of Significance"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Homogeneity
p5 <- ggplot(df_homo, aes(x = case, y = perc, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ comp, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Test of Homogeneity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Nullity
p6 <- ggplot(df_null, aes(x = case, y = perc, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ comp, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Test of Nullity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom")


# Figure S.3:
# Combine the plots into a 2x3 grid
pdf("figureS3.pdf", width = 13, height = 9)
(p1 | p2 | p3) / (p4 | p5 | p6) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()
