rm(list=ls())
library("mvtnorm")
library("glmnet")
library("Matrix")
library("multivar")
library("superheat")

# Source code:
source("multiVARSE_simulator.R")
source("multiVARSE_optimizer.R")
source("multiVARSE_debiasor.R")
source("multiVARSE_identifier.R")
source("multiVARSE_cv.R")
source("multiVARSE_tester.R")
source("multiVARSE_benchmark.R")
#################################################################################

K <- 10; d <- 20; p <- 1;
TT_vec <- rep(100,K)
list_data <- multiVARSE_simulator(K,d,p,TT_vec,s0=0.05,sk=0.05,eps=1)

# K <- 10; d <- 10; p <- 1;
# TT_vec <- rep(200,K)
# list_data <- multiVARSE_simulator(K,d,p,TT_vec,s0=0.05,sk=0.05,eps=1)

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
list_identified <- multiVARSE_identifier(list_debiased,eta_hat,thres_0_hat,K,d,p,TT_vec)

# Hypotheses tests on 
test_sig <- multiVARSE_tester(list_data,list_debiased,type="significance",K,d,p,TT_vec)
test_hetero <- multiVARSE_tester(list_data,list_debiased,type="heterogeneity",K,d,p,TT_vec)



# Benchmark method: multivar
list_benchmark <- multiVARSE_benchmark(list_data,K,d)


# Plots
par(mfrow=c(4,3))
matplot(cbind(list_identified$alpha_0_hat,list_data$alpha_0),
        type='l',ylab="",main="common,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[1]],list_data$alpha_k[[1]]),
        type='l',ylab="",main="unique1,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[2]],list_data$alpha_k[[2]]),
        type='l',ylab="",main="unique2,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[3]],list_data$alpha_k[[3]]),
        type='l',ylab="",main="unique3,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[4]],list_data$alpha_k[[4]]),
        type='l',ylab="",main="unique4,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[5]],list_data$alpha_k[[5]]),
        type='l',ylab="",main="unique5,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[6]],list_data$alpha_k[[6]]),
        type='l',ylab="",main="unique6,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[7]],list_data$alpha_k[[7]]),
        type='l',ylab="",main="unique7,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[8]],list_data$alpha_k[[8]]),
        type='l',ylab="",main="unique8,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[9]],list_data$alpha_k[[9]]),
        type='l',ylab="",main="unique9,multiVARSE")
matplot(cbind(list_identified$alpha_k_hat[[10]],list_data$alpha_k[[10]]),
        type='l',ylab="",main="unique10,multiVARSE")

round(c(norm(list_identified$alpha_0_hat-list_data$alpha_0,"2")/norm(list_data$alpha_0,"2"),
        norm(list_identified$alpha_k_hat[[1]] - list_data$alpha_k[[1]],"2")/norm(list_data$alpha_k[[1]],"2"),
        norm(list_identified$alpha_k_hat[[2]] - list_data$alpha_k[[2]],"2")/norm(list_data$alpha_k[[2]],"2"),
        norm(list_identified$alpha_k_hat[[3]] - list_data$alpha_k[[3]],"2")/norm(list_data$alpha_k[[3]],"2"),
        norm(list_identified$alpha_k_hat[[4]] - list_data$alpha_k[[4]],"2")/norm(list_data$alpha_k[[4]],"2"),
        norm(list_identified$alpha_k_hat[[5]] - list_data$alpha_k[[5]],"2")/norm(list_data$alpha_k[[5]],"2"),
        norm(list_identified$alpha_k_hat[[6]] - list_data$alpha_k[[6]],"2")/norm(list_data$alpha_k[[6]],"2"),
        norm(list_identified$alpha_k_hat[[7]] - list_data$alpha_k[[7]],"2")/norm(list_data$alpha_k[[7]],"2"),
        norm(list_identified$alpha_k_hat[[8]] - list_data$alpha_k[[8]],"2")/norm(list_data$alpha_k[[8]],"2"),
        norm(list_identified$alpha_k_hat[[9]] - list_data$alpha_k[[9]],"2")/norm(list_data$alpha_k[[9]],"2"),
        norm(list_identified$alpha_k_hat[[10]] - list_data$alpha_k[[10]],"2")/norm(list_data$alpha_k[[10]],"2")),3)


round(sum(list_data$alpha_0!=0 & list_identified$alpha_0_hat!=0)/sum(list_data$alpha_0!=0),3)
round(sum(list_data$alpha_0==0 & list_identified$alpha_0_hat==0)/sum(list_data$alpha_0==0),3)



# Plots, benchmark
par(mfrow=c(4,3))
matplot(cbind(as.vector(t(list_benchmark$common)),list_data$alpha_0),
        type='l',ylab="",main="common,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[1]])),list_data$alpha_k[[1]]),
        type='l',ylab="",main="unique1,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[2]])),list_data$alpha_k[[2]]),
        type='l',ylab="",main="unique2,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[3]])),list_data$alpha_k[[3]]),
        type='l',ylab="",main="unique3,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[4]])),list_data$alpha_k[[4]]),
        type='l',ylab="",main="unique4,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[5]])),list_data$alpha_k[[5]]),
        type='l',ylab="",main="unique5,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[6]])),list_data$alpha_k[[6]]),
        type='l',ylab="",main="unique6,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[7]])),list_data$alpha_k[[7]]),
        type='l',ylab="",main="unique7,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[8]])),list_data$alpha_k[[8]]),
        type='l',ylab="",main="unique8,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[9]])),list_data$alpha_k[[9]]),
        type='l',ylab="",main="unique9,multiVAR")
matplot(cbind(as.vector(t(list_benchmark$unique[[10]])),list_data$alpha_k[[10]]),
        type='l',ylab="",main="unique10,multiVAR")



round(c(norm(as.vector(t(list_benchmark$common))-list_data$alpha_0,"2")/norm(list_data$alpha_0,"2"),
        norm(as.vector(t(list_benchmark$unique[[1]])) - list_data$alpha_k[[1]],"2")/norm(list_data$alpha_k[[1]],"2"),
        norm(as.vector(t(list_benchmark$unique[[2]])) - list_data$alpha_k[[2]],"2")/norm(list_data$alpha_k[[2]],"2"),
        norm(as.vector(t(list_benchmark$unique[[3]])) - list_data$alpha_k[[3]],"2")/norm(list_data$alpha_k[[3]],"2"),
        norm(as.vector(t(list_benchmark$unique[[4]])) - list_data$alpha_k[[4]],"2")/norm(list_data$alpha_k[[4]],"2"),
        norm(as.vector(t(list_benchmark$unique[[5]])) - list_data$alpha_k[[5]],"2")/norm(list_data$alpha_k[[5]],"2"),
        norm(as.vector(t(list_benchmark$unique[[6]])) - list_data$alpha_k[[6]],"2")/norm(list_data$alpha_k[[6]],"2"),
        norm(as.vector(t(list_benchmark$unique[[7]])) - list_data$alpha_k[[7]],"2")/norm(list_data$alpha_k[[7]],"2"),
        norm(as.vector(t(list_benchmark$unique[[8]])) - list_data$alpha_k[[8]],"2")/norm(list_data$alpha_k[[8]],"2"),
        norm(as.vector(t(list_benchmark$unique[[9]])) - list_data$alpha_k[[9]],"2")/norm(list_data$alpha_k[[9]],"2"),
        norm(as.vector(t(list_benchmark$unique[[10]])) - list_data$alpha_k[[10]],"2")/norm(list_data$alpha_k[[10]],"2")),3)


round(sum(list_data$alpha_0!=0 & as.vector(t(list_benchmark$common))!=0)/sum(list_data$alpha_0!=0),3)
round(sum(list_data$alpha_0==0 & as.vector(t(list_benchmark$common))==0)/sum(list_data$alpha_0==0),3)
