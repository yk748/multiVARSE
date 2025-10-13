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
