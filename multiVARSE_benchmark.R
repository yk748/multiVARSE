###################################################################################
# The following is the bencmark methods for multiVAR models used in Section 5.
# The explanation of the method is included in Section 2.1. 
# The package loading (multivar) is included in the main loop of the simulator.
###################################################################################
multiVARSE_benchmark <- function(data_multivar,K,d,type){
  
  # ---------------------------------------------- #
  if (type == "adaptive"){
    multivar_model <- multivar::constructModel(data=data_multivar$data)
  }else if (type == "non-adaptive"){
    multivar_model <- multivar::constructModel(data=data_multivar$data,lassotype="standard")
  }
  
  cv_multivar <- multivar::cv.multivar(multivar_model)
  return(cv_multivar$mats)
}
