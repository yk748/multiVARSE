
get_sim <- function(p,K,d,TT_mean,case,iter_max=50){
  
  set.seed(123)
  if (TT_mean == 50){
    TT <- sample(seq(45,55),K,replace=TRUE)
  }else if (TT_mean == 200){
    TT <- sample(seq(190,210),K,replace=TRUE)
  }
  
  if (case == "case1"){
    s0 <- 0.02; sk <- rep(0.04,K)
  }else if (case == "case2") {
    s0 <- 0.03; sk <- rep(0.03,K)
  }else if (case == "case3"){
    s0 <- 0.04; sk <- rep(0.02,K)
  }
  
  
  library("doParallel")
  library("dqrng")
  library("doRNG")
  
  list_dummy <- list()
  cl <- makeCluster(50)
  registerDoParallel(cl)
  
  sim_multivar <- foreach(iter = 1:iter_max,
                        .errorhandling = 'remove',
                        .packages = c("mvtnorm","glmnet","Matrix")) %dopar% {
                          
                          source("multiVARSE_simulator.R")
                          source("multiVARSE_optimizer.R")
                          source("multiVARSE_debiasor.R")
                          source("multiVARSE_cv.R")
                          source("multiVARSE_identifier.R")
                          source("multiVARSE_tester.R")
                          source("multiVARSE_test_evaluator.R")
                          source("multiVARSE_benchmark.R")
                          
                          #####################################################################
                          # Data generation:
                          list_data <- multiVARSE_simulator(K,d,p,TT,s0=s0,sk=sk,eps=1)
                          
                          #####################################################################
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
                          
                          ###################################################################
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
                          
                          #####################################################################
                          # Estimation:

                          start_method <- proc.time()
                          list_opt <- list()
                          list_debiased <- list()
                          list_identified <- list()
                          for (k in 1:K){
                            # Individual estimation:
                            X_t <- list_data$x_t[[k]]
                            list_opt[[k]] <- multiVARSE_optimizer(X_t,d,p,TT[k])

                            # Individual debiasing:
                            beta_hat <- list_opt[[k]]$beta_hat
                            list_debiased[[k]] <- multiVARSE_debiasor(beta_hat,X_t,d,p,TT[k])
                            cat(k,"th subject done\n")
                          }

                          # Cross validation for thresholds
                          cv_tuning <- multiVARSE_cv2(list_data$x_t,list_debiased,list_opt,
                                                      K,d,p,TT,thres_type="hard",
                                                      nfold=5,grid_eta=10,grid_thres_0=10,grid_thres_k=5)

                          # Identifying common and unique components
                          eta_hat <- cv_tuning$eta_hat
                          thres_0_hat <- cv_tuning$thres_0_hat
                          thres_k_hat <- cv_tuning$thres_k_hat
                          list_identified <- multiVARSE_identifier2(list_debiased,eta_hat,thres_0_hat,
                                                                    thres_k_hat,thres_type="hard",
                                                                    K,d,p,TT)
                          time_method <- proc.time() - start_method
                          
                          ####################################################################
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
                          
                          #####################################################################
                          # Hypotheses tests:
                          test_null <- multiVARSE_tester(list_data$x_t,list_debiased,list_opt,list_identified,
                                                         type="nullity",alpha=0.05,K,d,p,TT)
                          test_homo <- multiVARSE_tester(list_data$x_t,list_debiased,list_opt,list_identified,
                                                         type="homogeneity",alpha=0.05,K,d,p,TT)
                          test_sig <- multiVARSE_tester(list_data$x_t,list_debiased,list_opt,list_identified,
                                                        type="significance",alpha=0.05,K,d,p,TT)
                          
                          #####################################################################
                          # Evaluation - Hypothesis tests:
                          eval_null <- multiVARSE_test_evaluator(list_data,test_null,K,d,p,
                                                                 type="nullity",alpha=0.05)
                          eval_homo <- multiVARSE_test_evaluator(list_data,test_homo,K,d,p,
                                                                 type="homogeneity",alpha=0.05)
                          eval_sig <- multiVARSE_test_evaluator(list_data,test_sig,K,d,p,
                                                                type="significance",alpha=0.05)

                          ######################################################################
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
                          output$eval_test <- c(eval_null$FDR,
                                                eval_homo$FDR,
                                                eval_sig$FDR,
                                                eval_null$power,
                                                eval_homo$power,
                                                eval_sig$power)
                          list_dummy[[iter]] <- output
                        }  
  stopCluster(cl)
  print(length(sim_multivar))
  if (length(sim_multivar) == 0) {
    warning("No successful iterations were returned.")
  }
  save(sim_multivar,
       file=paste0("multivar_d",d,"_T",TT_mean,"_subject",K,"_",case,".RData"))
}

#####################################################################################
p <- 1
K <- 10
# K <- 15

d_list <- c(10,20)
TT_list <- c(50,200)
K_list <- c(10,15)
case_list <- c("case1","case2","case3")

for (i in 1:3){
  d <- d_list[i]
  
  for (j in 1:2){
    TT_mean <- TT_list[j]
      
    for (k in 1:2){
      K <- K_list[k]
      
      for (l in 2:3){
        case <- case_list[l]
        
        cat("The current iteration is d:",d,",TT:",TT_mean,",subject:",K,",case:",case,"\n")
        get_sim(p,K,d,TT_mean,case,iter_max=50)
        cat("Done.\n")
      } 
    }
  }
}


###################################################################################
library(ggplot2)
library(latex2exp)
library(patchwork)
library(Matrix)
library(RColorBrewer)
 
#####################################################################################
# Estimation:
p <- 1

d_list <- c(10,20)
TT_list <- c(50,200)
K_list <- c(10,15)
case_list <- c("case1","case2","case3")
iter_max <- 50

df_est_prop <- data.frame()
df_est_ls <- data.frame()
df_est_als <- data.frame()

for (i in 1:2){
  d <- d_list[i]

  for (j in 1:2){
    TT_mean <- TT_list[j]

    for (k in 1:2){
      K <- K_list[k]

      for (l in 1:3){
        case <- case_list[l]

        cat("The current iteration is d:",d,",TT:",TT_mean,",subject:",K,",case:",case,"\n")
        load(paste0("multivar_d",d,"_T",TT_mean,"_subject",K,"_",case,".RData"))

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
                                   K = rep(K,iter_max),
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
                                      K = rep(K,iter_max),
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
                                      K = rep(K,iter_max),
                                      case = rep(case,iter_max),
                                      method = rep("als",iter_max)))
      }
    }
  }
}



# --------------------------------------------------------------- #
df_est <- rbind(df_est_prop,df_est_ls,df_est_als)
df_est$d <- factor(df_est$d,levels=unique(df_est$d),
                   labels=c("d=10","d=20"))
df_est$TT <- factor(df_est$TT,levels=unique(df_est$TT),
                    labels=c("T=50","T=200"))
df_est$K <- factor(df_est$K,levels=unique(df_est$K),
                   labels=c("K=10","K=15"))
df_est$case <- factor(df_est$case,levels=unique(df_est$case),
                       labels=c("High","Med","Low"))
df_est$method <- factor(df_est$method,levels=unique(df_est$method),
                        labels=c("Proposed","Multi-VAR","Multi-VAR (A)"))


df_rmse_comm <- df_est[,c(1,10,11,12,13,14)]
df_rmse_uniq <- df_est[,c(2,10,11,12,13,14)]
df_rmse_tot <- df_est[,c(3,10,11,12,13,14)]
df_sens_comm <- df_est[,c(4,10,11,12,13,14)]
df_sens_uniq <- df_est[,c(5,10,11,12,13,14)]
df_sens_tot <- df_est[,c(6,10,11,12,13,14)]
df_spec_comm <- df_est[,c(7,10,11,12,13,14)]
df_spec_uniq <- df_est[,c(8,10,11,12,13,14)]
df_spec_tot <- df_est[,c(9,10,11,12,13,14)]


# RMSE:
p1 <- ggplot(df_rmse_comm, aes(x = case, y = rmse_comm, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("RMSE of $\\alpha^{(0)}$"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

p2 <- ggplot(df_rmse_uniq, aes(x = case, y = rmse_uniq, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Mean RMSE of $\\alpha^{(k)}$"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom")

# Sens:
p3 <- ggplot(df_sens_comm, aes(x = case, y = sens_comm, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Sens of $\\alpha^{(0)}$"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

p4 <- ggplot(df_sens_uniq, aes(x = case, y = sens_uniq, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Mean Sens of $\\alpha^{(k)}$"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom")


# Spec:
p5 <- ggplot(df_spec_comm, aes(x = case, y = spec_comm, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Spec of $\\alpha^{(0)}$"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

p6 <- ggplot(df_spec_uniq, aes(x = case, y = spec_uniq, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Mean Spec of $\\alpha^{(k)}$"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set1", name = "") +
  theme(legend.position = "bottom")



# Combine the plots into a 3x2 grid
pdf("est.pdf", width = 13, height = 13)
(p1 | p2) / (p3 | p4) / (p5 | p6) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()


#####################################################################################
# Time:
p <- 1

d_list <- c(10,20)
TT_list <- c(50,200)
K_list <- c(10,15)
case_list <- c("case1","case2","case3")
iter_max <- 50

df_time <- data.frame()

for (i in 1:2){
  d <- d_list[i]

  for (j in 1:2){
    TT_mean <- TT_list[j]

    for (k in 1:2){
      K <- K_list[k]

      for (l in 1:3){
        case <- case_list[l]

        cat("The current iteration is d:",d,",TT:",TT_mean,",subject:",K,",case:",case,"\n")
        load(paste0("multivar_d",d,"_T",TT_mean,"_subject",K,"_",case,".RData"))

        # --------------------------------------------- #
        df_time_tmp <- c(mapply(iter=1:iter_max,function(iter)
          sim_multivar[[iter]]$time$ls[3]/sim_multivar[[iter]]$time$method[3]),
                         mapply(iter=1:iter_max,function(iter)
          sim_multivar[[iter]]$time$als[3]/sim_multivar[[iter]]$time$method[3]))
        df_time <- rbind(df_time,
                         data.frame(value=df_time_tmp,
                                    d = rep(d,(2*iter_max)),
                                    TT = rep(TT_mean,(2*iter_max)),
                                    K = rep(K,(2*iter_max)),
                                    case = rep(case,(2*iter_max)),
                                    method = c(rep("ls",iter_max),rep("als",iter_max))))
      }
    }
  }
}

# --------------------------------------------------------------- #
df_time$d <- factor(df_time$d,levels=unique(df_time$d),
                   labels=c("d=10","d=20"))
df_time$TT <- factor(df_time$TT,levels=unique(df_time$TT),
                    labels=c("T=50","T=200"))
df_time$K <- factor(df_time$K,levels=unique(df_time$K),
                   labels=c("K=10","K=15"))
df_time$case <- factor(df_time$case,levels=unique(df_time$case),
                      labels=c("High","Med","Low"))
df_time$method <- factor(df_time$method,levels=unique(df_time$method),
                        labels=c("Multi-VAR","Multi-VAR (A)"))

pdf("cputime.pdf", width = 13, height = 13)
ggplot(df_time, aes(x = case, y = value, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = "CPU Time Ratio (Benchmark / Proposed)",
       x = "Heterogeneity levels") +
  scale_fill_manual(name = "", values=brewer.pal(3, "Set1")[2:3]) +
  theme(legend.position = "bottom")
dev.off()

#####################################################################################
# Hypothesis tests:
p <- 1

d_list <- c(10,20)
TT_list <- c(50,200)
K_list <- c(10,15)
case_list <- c("case1","case2","case3")
iter_max <- 50

df_test <- data.frame()

for (i in 1:2){
  d <- d_list[i]

  for (j in 1:2){
    TT_mean <- TT_list[j]

    for (k in 1:2){
      K <- K_list[k]

      for (l in 1:3){
        case <- case_list[l]

        cat("The current iteration is d:",d,",TT:",TT_mean,",subject:",K,",case:",case,"\n")
        load(paste0("multivar_d",d,"_T",TT_mean,"_subject",K,"_",case,".RData"))

        df_test_tmp <- data.frame(t(mapply(iter=1:iter_max,function(iter)
          sim_multivar[[iter]]$eval_test)))
        names(df_test_tmp) <- c("FDR_null","FDR_homo","FDR_sig",
                                "Power_null","Power_homo","Power_sig")
        df_test <- rbind(df_test,
                         data.frame(df_test_tmp,
                                    d = rep(d,iter_max),
                                    TT = rep(TT_mean,iter_max),
                                    K = rep(K,iter_max),
                                    case = rep(case,iter_max)))

      }
    }
  }
}
df_test <- na.omit(df_test)

# --------------------------------------------------------------- #
df_test$d <- factor(df_test$d,levels=unique(df_test$d),
                    labels=c("d=10","d=20"))
df_test$TT <- factor(df_test$TT,levels=unique(df_test$TT),
                     labels=c("T=50","T=200"))
df_test$K <- factor(df_test$K,levels=unique(df_test$K),
                    labels=c("K=10","K=15"))
df_test$case <- factor(df_test$case,levels=unique(df_test$case),
                       labels=c("High","Med","Low"))



df_fdr_null <- df_test[,c(1,7,8,9,10)]
df_fdr_homo <- df_test[,c(2,7,8,9,10)]
df_fdr_sig <- df_test[,c(3,7,8,9,10)]
df_power_null <- df_test[,c(4,7,8,9,10)]
df_power_homo <- df_test[,c(5,7,8,9,10)]
df_power_sig <- df_test[,c(6,7,8,9,10)]


# Plot 7
p7 <- ggplot(df_fdr_null, aes(x = case, y = FDR_null, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("FDR of Test of Nullity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Plot 8
p8 <- ggplot(df_power_null, aes(x = case, y = Power_null, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Power of Test of Nullity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom")

# Plot 9
p9 <- ggplot(df_fdr_homo, aes(x = case, y = FDR_homo, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("FDR of Test of Homogeneity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Plot 10
p10 <- ggplot(df_power_homo, aes(x = case, y = Power_homo, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Power of Test of Homogeneity"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom")

# Plot 11
p11 <- ggplot(df_fdr_sig, aes(x = case, y = FDR_sig, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("FDR of Test of Significance"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom",
        strip.text.y = element_blank(),
        strip.background = element_blank())

# Plot 12
p12 <- ggplot(df_power_sig, aes(x = case, y = Power_sig, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.85)) +
  facet_grid(d+TT ~ K, scales = "free") +
  theme_minimal() +
  labs(y = TeX("Power of Test of Significance"),
       x = "Heterogeneity levels") +
  scale_fill_brewer(palette = "Set2", name = "") +
  theme(legend.position = "bottom")

# Combine the plots into a 3x2 grid
pdf("test.pdf", width = 13, height = 13)
(p7 | p8) / (p9 | p10) / (p11 | p12) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()

