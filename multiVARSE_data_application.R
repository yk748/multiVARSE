###################################################################################
# The following is the data application procedure in Section 6.
# The final result is diplayed in Figure 3 in Section 6.2.
###################################################################################
rm(list=ls())

library("mvtnorm")
library("glmnet")
library("Matrix")
library("multivar")
library("magic")
library("ggplot2")
library("reshape2")
library("patchwork")

# Source code:
source("multiVARSE_optimizer.R")
source("multiVARSE_debiasor.R")
source("multiVARSE_cv.R")
source("multiVARSE_identifier.R")
source("multiVARSE_benchmark.R")
source("multiVARSE_tester.R")

###################################################################################
# The following is for data configuration & preprocessing
# The corresponding codes are commentized as the processed data is available 
# on the author's GitHub.
###################################################################################
# # --------------------------------------------------------------- #
# # Workbench directory:
# library(ciftiTools)
# ciftiTools.setOption("wb_path", "C:/workbench/bin_windows64/wb_command.exe")
# 
# subj <- c(100307, 117122, 162329, 182739, 192540,
#           194140, 199150, 255639, 499566, 729557,
#           826353, 901139)
# 
# yeo17_names <- c(
#   "Visual_A", "Visual_B", "Somatomotor_A", "Somatomotor_B", 
#   "DorsalAttention_A", "DorsalAttention_B", "VentralAttention_A", "VentralAttention_B", 
#   "Limbic_A", "Limbic_B", "Frontoparietal_A", "Frontoparietal_B", 
#   "Frontoparietal_C", "DefaultMode_A", "DefaultMode_B", "DefaultMode_C", "DefaultMode_D"
# )
# 
# # --------------------------------------------------------------- #
# # Prepare Atlas & Network Mapping
# atlas_file  <- "Schaefer2018_400Parcels_17Networks_order.dlabel.nii"   
# atlas_cifti <- read_cifti(atlas_file)
# 
# # Safe vector conversion using as.matrix()
# atlas_labels_vec <- as.numeric(as.matrix(atlas_cifti))
# 
# # Robust Manual Mapping for Schaefer 400 (17 Networks)
# # Official number of parcels assigned to each of the 17 networks per hemisphere:
# parcels_per_net_LH <- c(9, 22, 14, 20, 12, 17, 10, 12, 4, 4, 11, 14, 5, 14, 19, 7, 6)
# parcels_per_net_RH <- c(11, 21, 11, 23, 13, 14, 11, 10, 4, 5, 11, 14, 5, 16, 21, 6, 5)
# 
# lh_mapping <- rep(yeo17_names, times = parcels_per_net_LH)
# rh_mapping <- rep(yeo17_names, times = parcels_per_net_RH)
# 
# # Combine Left and Right hemispheres sequentially (Total 400 strings)
# parcel_network_mapping <- c(lh_mapping, rh_mapping)
# 
# # --------------------------------------------------------------- #
# # Extraction Function
# extract_roi_ts_cifti <- function(fmri_matrix, atlas_labels_vec) {
#   # Unique spatial parcel IDs present in the file (ignoring background/0)
#   roi_ids <- unique(atlas_labels_vec)
#   roi_ids <- roi_ids[roi_ids != 0 & !is.na(roi_ids)]
#   roi_ids <- sort(roi_ids) # Ensure order aligns to 1:400
#   
#   nTime <- ncol(fmri_matrix) 
#   nROIs  <- length(roi_ids)
#   
#   roi_ts <- matrix(NA, nrow = nTime, ncol = nROIs)
#   
#   for (i in seq_along(roi_ids)) {
#     roi <- roi_ids[i]
#     idx <- which(atlas_labels_vec == roi)
#     if(length(idx) == 0) next
#     
#     # Rows are space, Columns are time in ciftiTools matrix output
#     roi_ts[, i] <- colMeans(fmri_matrix[idx, , drop = FALSE], na.rm = TRUE)
#   }
#   
#   colnames(roi_ts) <- paste0("ROI_", roi_ids)
#   return(roi_ts)
# }
# 
# # --------------------------------------------------------------- #
# # Main Subject Loop (Fixed with Guardrails)
# for (i in 1:length(subj)){
#   
#   # File paths
#   phase1_file <- paste0(subj[i],"_3T_tfMRI_EMOTION_preproc/",subj[i],"/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_Atlas_MSMAll.dtseries.nii")
#   phase2_file <- paste0(subj[i],"_3T_tfMRI_EMOTION_preproc/",subj[i],"/MNINonLinear/Results/tfMRI_EMOTION_RL/tfMRI_EMOTION_RL_Atlas_MSMAll.dtseries.nii")
#   
#   # Skip if files don't exist
#   if(!file.exists(phase1_file) || !file.exists(phase2_file)) {
#     cat("Missing files for subject", subj[i], "- skipping.\n")
#     next
#   }
#   
#   # Load and cleanly matrix-convert fMRI data
#   phase1_matrix <- as.matrix(read_cifti(phase1_file))
#   phase2_matrix <- as.matrix(read_cifti(phase2_file))  
#   
#   # Extract individual parcel time series (Time x ROIs Matrix)
#   roi_ts_phase1 <- extract_roi_ts_cifti(phase1_matrix, atlas_labels_vec)
#   roi_ts_phase2 <- extract_roi_ts_cifti(phase2_matrix, atlas_labels_vec)
#   
#   # Concatenate LR and RL runs chronologically 
#   roi_ts_combined <- rbind(roi_ts_phase1, roi_ts_phase2) 
#   
#   # Aggregate the Parcels into the 17 Yeo Networks
#   network_ts <- matrix(NA, nrow = nrow(roi_ts_combined), ncol = length(yeo17_names))
#   colnames(network_ts) <- yeo17_names
#   
#   for (n in 1:length(yeo17_names)) {
#     net <- yeo17_names[n]
#     
#     # 1. Find which generic parcel IDs (1 to 400) belong to this network string
#     matching_parcel_ids <- which(parcel_network_mapping == net)
#     
#     # 2. Convert those IDs to the actual column names expected ("ROI_1", "ROI_2", etc.)
#     expected_col_names <- paste0("ROI_", matching_parcel_ids)
#     
#     # 3. FIX: Safely match against the columns that ACTUALLY exist in the data matrix
#     col_idx <- match(expected_col_names, colnames(roi_ts_combined))
#     col_idx <- col_idx[!is.na(col_idx)] # Remove any missing/out-of-bounds indices
#     
#     # 4. Calculate mean only if columns were safely found
#     if(length(col_idx) > 0) {
#       network_ts[, net] <- rowMeans(roi_ts_combined[, col_idx, drop = FALSE], na.rm = TRUE)
#     } else {
#       warning(paste("No matching data columns found for network:", net, "in subject", subj[i]))
#     }
#   }
#   
#   # Save Network time series
#   write.csv(network_ts, file=paste0("new_HCP_Emotion_",subj[i],".csv"), row.names=FALSE)
#   cat(i, "th subject (", subj[i], ") done.\n", sep="")
# }


###################################################################################
# The following is where the main part begins.
# First, we load the processed data.
###################################################################################
yeo17_names <- c(
  "Visual_A",      # 1
  "Visual_B",      # 2
  "Somatomotor_A", # 3
  "Somatomotor_B", # 4
  "DorsalAttention_A", # 5
  "DorsalAttention_B", # 6
  "VentralAttention_A", # 7
  "VentralAttention_B", # 8
  "Limbic_A",      # 9
  "Limbic_B",      # 10
  "Frontoparietal_A", # 11
  "Frontoparietal_B", # 12
  "Frontoparietal_C", # 13
  "DefaultMode_A", # 14
  "DefaultMode_B", # 15
  "DefaultMode_C", # 16
  "DefaultMode_D"  # 17
)

K_idx <- c(100307,117122,162329,182739,192540,
           194140,199150,255639,499566,729557,
           826353,901139)
K <- length(K_idx)
d <- 17
p <- 1
TT <- rep(165,K)

list_data <- list()
for (k in 1:K){
  
  ts <- read.csv(paste0("./HCP_Emotion_",K_idx[k],".csv"),header=TRUE)
  x_t <- scale(ts)
  list_data$x_t[[k]] <- t(x_t[-c(c(1:5),c(171:176)),])
}


###################################################################################
# We apply the benchmark methods.
###################################################################################
# Benchmark method: multivar
data_multivar <- list()
for (k in 1:K){
  data_multivar[["data"]][[k]] <- t(list_data$x_t[[k]])
}

list_benchmark_ls <- multiVARSE_benchmark(data_multivar,K,d,type="non-adaptive")
list_benchmark_als <- multiVARSE_benchmark(data_multivar,K,d,type="adaptive")


###################################################################################
# We apply the proposed estimation method.
###################################################################################
# Estimation:
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

# Identification
# Cross validation for thresholds
cv_tuning <- multiVARSE_cv(list_data$x_t,list_debiased,list_opt,
                            agg_type = "redescending",
                            K,d,p,TT,thres_type="hard",nfold=5)

# Identifying common and unique components
eta_hat <- cv_tuning$eta_hat
thres_0_hat <- cv_tuning$thres_0_hat
thres_k_hat <- cv_tuning$thres_k_hat
list_identified <- multiVARSE_identifier(list_debiased,
                                          agg_type = "redescending",
                                          eta_hat,thres_0_hat,thres_k_hat,
                                          thres_type="hard",K,d,p,TT)

###################################################################################
# We apply the proposed testing method.
###################################################################################
# Hypotheses tests:
test_sig <- multiVARSE_tester(list_data$x_t,list_debiased,
                              list_opt,list_identified,
                              type="significance",alpha=0.05,K,d,p,TT)
test_homo <- multiVARSE_tester(list_data$x_t,list_debiased,
                               list_opt,list_identified,
                               type="homogeneity",alpha=0.05,K,d,p,TT)
test_null <- multiVARSE_tester(list_data$x_t,list_debiased,
                               list_opt,list_identified,
                               type="nullity",alpha=0.05,K,d,p,TT)


###################################################################################
# We begin processing the results for display.
###################################################################################
list_result <- list()

# Estimation:
list_result$mat_com <- t(matrix(list_identified$alpha_0_hat,d,d,byrow=FALSE))
list_result$mat_uniq <- matrix(0,d,d)
list_result$mat_ind <- matrix(0,d,d)

for (k in 1:K){
  list_result$mat_uniq <- list_result$mat_uniq + 1*(t(matrix(list_identified$alpha_k_hat[[k]],d,d,byrow=FALSE))!=0)
  list_result$mat_ind <- list_result$mat_ind + 1*(t(matrix(list_identified$alpha_k_hat[[k]]+list_result$mat_com,d,d,byrow=FALSE))!=0)
}
colnames(list_result$mat_com) <- rownames(list_result$mat_com) <- yeo17_names
colnames(list_result$mat_uniq) <- rownames(list_result$mat_uniq) <- yeo17_names
colnames(list_result$mat_ind) <- rownames(list_result$mat_ind) <- yeo17_names


# Testing:
list_result$mat_com_filt <- list_result$mat_com* (1*t(matrix(test_sig$decision,d,d)))
list_result$mat_uniq_filt <- list_result$mat_uniq* (1*t(matrix(test_homo$decision,d,d)))
list_result$mat_ind_filt <- list_result$mat_ind* (1*t(matrix(test_null$decision,d,d)))

colnames(list_result$mat_com_filt) <- rownames(list_result$mat_com_filt) <- yeo17_names
colnames(list_result$mat_uniq_filt) <- rownames(list_result$mat_uniq_filt) <- yeo17_names
colnames(list_result$mat_ind_filt) <- rownames(list_result$mat_ind_filt) <- yeo17_names


# BM - Lasso:
list_result$mat_com_ls <- t(matrix(as.vector(t(list_benchmark_ls$common)),d,d,byrow=FALSE))
list_result$mat_uniq_ls <- matrix(0,d,d)
list_result$mat_ind_ls <- matrix(0,d,d)
for (k in 1:K){
  list_result$mat_uniq_ls <- list_result$mat_uniq_ls + 1*(t(matrix(list_benchmark_ls$unique[[k]],d,d,byrow=FALSE))!=0)
  list_result$mat_ind_ls <- list_result$mat_ind_ls + 1*(t(matrix(list_benchmark_ls$unique[[k]]+list_result$mat_com_ls,d,d,byrow=FALSE))!=0)
}
colnames(list_result$mat_com_ls) <- rownames(list_result$mat_com_ls) <- yeo17_names
colnames(list_result$mat_uniq_ls) <- rownames(list_result$mat_uniq_ls) <- yeo17_names
colnames(list_result$mat_ind_ls) <- rownames(list_result$mat_ind_ls) <- yeo17_names

# BM - Adaptive Lasso:
list_result$mat_com_als <- t(matrix(as.vector(t(list_benchmark_als$common)),d,d,byrow=FALSE))
list_result$mat_uniq_als <- matrix(0,d,d)
list_result$mat_ind_als <- matrix(0,d,d)
for (k in 1:K){
  list_result$mat_uniq_als <- list_result$mat_uniq_als + 1*(t(matrix(list_benchmark_als$unique[[k]],d,d,byrow=FALSE))!=0)
  list_result$mat_ind_als <- list_result$mat_ind_als + 1*(t(matrix(list_benchmark_als$unique[[k]]+list_result$mat_com_als,d,d,byrow=FALSE))!=0)
}
colnames(list_result$mat_com_als) <- rownames(list_result$mat_com_als) <- yeo17_names
colnames(list_result$mat_uniq_als) <- rownames(list_result$mat_uniq_als) <- yeo17_names
colnames(list_result$mat_ind_als) <- rownames(list_result$mat_ind_als) <- yeo17_names


# --------------------------------------------------------------- #
# Common paths:
round(sum(list_result$mat_com_ls!=0 & list_result$mat_com!=0)/sum(list_result$mat_com!=0),3)
round(sum(list_result$mat_com_als!=0 & list_result$mat_com!=0)/sum(list_result$mat_com!=0),3)

round(sum(list_result$mat_com_ls!=0 & list_result$mat_com_filt!=0)/sum(list_result$mat_com_filt!=0),3)
round(sum(list_result$mat_com_als!=0 & list_result$mat_com_filt!=0)/sum(list_result$mat_com_filt!=0),3)

# Unique paths (frequency):
round(sum(list_result$mat_uniq_ls!=0)/d^2,3)
round(sum(list_result$mat_uniq_als!=0)/d^2,3)
round(sum(list_result$mat_uniq!=0)/d^2,3)
round(sum(list_result$mat_uniq_filt!=0)/d^2,3)

summary(as.vector(list_result$mat_uniq_ls))
summary(as.vector(list_result$mat_uniq_als))
summary(as.vector(list_result$mat_uniq))
summary(as.vector(list_result$mat_uniq))


# Individual paths (frequency):
round(sum(list_result$mat_ind_ls!=0)/d^2,3)
round(sum(list_result$mat_ind_als!=0)/d^2,3)
round(sum(list_result$mat_ind!=0)/d^2,3)
round(sum(list_result$mat_ind_filt!=0)/d^2,3)

summary(as.vector(list_result$mat_ind_ls))
summary(as.vector(list_result$mat_ind_als))
summary(as.vector(list_result$mat_ind))
summary(as.vector(list_result$mat_ind_filt))

# --------------------------------------------------------------- #
# BM - Lasso:
df_common_ls <- {
  mat <- list_result$mat_com_ls
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_unique_ls <- {
  mat <- list_result$mat_uniq_ls
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_ind_ls <- {
  mat <- list_result$mat_ind_ls
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_common_ls$rows <- factor(df_common_ls$rows, levels = rev(yeo17_names))
df_unique_ls$rows <- factor(df_unique_ls$rows, levels = rev(yeo17_names))
df_ind_ls$rows <- factor(df_ind_ls$rows, levels = rev(yeo17_names))

# BM - Adaptive Lasso:
df_common_als <- {
  mat <- list_result$mat_com_als
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_unique_als <- {
  mat <- list_result$mat_uniq_als
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_ind_als <- {
  mat <- list_result$mat_ind_als
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_common_als$rows <- factor(df_common_als$rows, levels = rev(yeo17_names))
df_unique_als$rows <- factor(df_unique_als$rows, levels = rev(yeo17_names))
df_ind_als$rows <- factor(df_ind_als$rows, levels = rev(yeo17_names))

# Estimation:
df_common <- {
  mat <- list_result$mat_com
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_unique <- {
  mat <- list_result$mat_uniq
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_ind <- {
  mat <- list_result$mat_ind
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_common$rows <- factor(df_common$rows, levels = rev(yeo17_names))
df_unique$rows <- factor(df_unique$rows, levels = rev(yeo17_names))
df_ind$rows <- factor(df_ind$rows, levels = rev(yeo17_names))

# Estimation + Hypothesis test:
df_common_filt <- {
  mat <- list_result$mat_com_filt
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_unique_filt <- {
  mat <- list_result$mat_uniq_filt
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_ind_filt <- {
  mat <- list_result$mat_ind_filt
  colnames(mat) <- rownames(mat) <- yeo17_names
  setNames(melt(mat), c("rows", "vars", "values"))
}
df_common_filt$rows <- factor(df_common_filt$rows, levels = rev(yeo17_names))
df_unique_filt$rows <- factor(df_unique_filt$rows, levels = rev(yeo17_names))
df_ind_filt$rows <- factor(df_ind_filt$rows, levels = rev(yeo17_names))



###################################################################################
# We begin creating Figure 3 in Section 6.2
###################################################################################
# common setting:
red <- rgb(255, 0, 90, maxColorValue = 255)
green <- rgb(90, 168, 0, maxColorValue = 255)
blue <- rgb(0, 152, 233, maxColorValue = 255)
yellow <- rgb(242, 147, 24, maxColorValue = 255)

fore <- "white"
text_color <- rgb(51, 51, 51, maxColorValue = 255)
grid_color <- rgb(51, 51, 51, maxColorValue = 255)
plot_background <- "white"


limit_val <- c(-1,1)
colval_low <- colorRampPalette(c(red, fore))
colval_high <- colorRampPalette(c(fore, blue))
colors_val <- c(colval_low(6)[1:3], fore, colval_high(6)[4:6])

limit_cnt <- c(0,K)
colcnt_low <- colorRampPalette(c(yellow, fore))
colcnt_high <- colorRampPalette(c(fore, green))
colors_cnt <- c(fore, colcnt_low(6)[1:3], colcnt_high(6)[4:6])


# Shared scales for fill
shared_fill_scale_val <- scale_fill_gradientn(
  name = "Values",
  colors = colors_val,
  limits = limit_val,
  na.value = fore,
  guide = guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black",
    ticks.linewidth = 0.1,
    frame.linewidth = 0.1
  )
)

shared_fill_scale_cnt <- scale_fill_gradientn(
  name = "Counts",
  colors = colors_cnt,
  limits = limit_cnt,
  na.value = fore,
  guide = guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black",
    ticks.linewidth = 0.1,
    frame.linewidth = 0.1
  )
)

# Shared tile plot theme
common_tile_theme <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = text_color),
  axis.text.y = element_text(size = 6, color = text_color),
  axis.ticks = element_blank(),
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
  strip.text = element_text(hjust = 0, color = text_color, size = 8, face = "bold"),
  strip.background = element_rect(fill = plot_background, color = plot_background),
  panel.spacing.x = unit(0.3, "cm"),
  panel.spacing.y = unit(0.3, "cm"),
  legend.background = element_rect(fill = plot_background, color = plot_background),
  legend.title = element_text(size = 8, color = text_color, face = "bold"),
  legend.title.align = 1,
  legend.text.align = 1,
  legend.text = element_text(size = 8, color = text_color),
  plot.background = element_rect(fill = plot_background, color = plot_background),
  plot.margin = margin(1, 1, 1, 1, "mm"),
  plot.title = element_text(hjust = 0, color = text_color, face = "bold", 
                            vjust = 1.2),  # hjust=0.5 centers title
  axis.title.x = element_text(vjust = -0.2),
  axis.title.y = element_text(vjust = 0.5)
)

no_y_axis <- theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank()
)


# --------------------------------------------------------------- #
# First row plots: Common paths
plot_common_ls <- ggplot(df_common_ls, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_val +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Common paths - multi-VAR") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)])

plot_common_als <- ggplot(df_common_als, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_val +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Common paths - multi-VAR(A)") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

plot_common <- ggplot(df_common, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_val +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Common paths - Estimation") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

plot_common_filt <- ggplot(df_common_filt, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_val +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Common paths - Tests") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

# Second row plots: Unique paths
plot_unique_ls <- ggplot(df_unique_ls, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Unique paths - multi-VAR") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) 

plot_unique_als <- ggplot(df_unique_als, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Unique paths - multi-VAR(A)") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

plot_unique <- ggplot(df_unique, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Unique paths - Estimation") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

plot_unique_filt <- ggplot(df_unique_filt, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "Frequency", x = NULL, y = NULL, title = "Unique paths - Tests") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

# Third row plots: Individual paths
plot_ind_ls <- ggplot(df_ind_ls, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Individual paths - multi-VAR") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) 

plot_ind_als <- ggplot(df_ind_als, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Individual paths - multi-VAR(A)") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

plot_ind <- ggplot(df_ind, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Individual paths - Estimation") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis

plot_ind_filt <- ggplot(df_ind_filt, aes(y = rows, x = vars, fill = values)) +
  geom_tile() +
  shared_fill_scale_cnt +
  coord_equal() +
  common_tile_theme +
  labs(fill = "", x = NULL, y = NULL, title = "Individual paths - Tests") +
  scale_x_discrete(breaks = yeo17_names[seq(1, length(yeo17_names), by = 2)]) +
  scale_y_discrete(breaks = yeo17_names[seq(2, length(yeo17_names), by = 2)]) +
  no_y_axis


# --------------------------------------------------------------- #
# First row (Common paths) - legend on the right
row_common <- (plot_common_ls | plot_common_als | plot_common | plot_common_filt) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Second row (Unique paths) - legend on the right
row_unique <- (plot_unique_ls | plot_unique_als | plot_unique | plot_unique_filt) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Third row (Individual paths) - legend on the right
row_ind <- (plot_ind_ls | plot_ind_als | plot_ind | plot_ind_filt) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Stack rows vertically
final_plot <- (row_common / row_unique / row_ind) +
  plot_layout(
    ncol = 1,
    heights = c(2.5, 2.5, 2.5)  # adjust row heights
  ) &
  theme(plot.title = element_text(size = 12))


# Figures 3:
pdf("figure3.pdf", width = 13, height = 13)
final_plot
dev.off()

