#' @title Model Selection
#' @description Model selection based on AIC and BIC information criteria for models in IMIX
#'
#' @param loglik Full log likelihood, result output from IMIX or a sub model in IMIX: `Full MaxLogLik final`
#' @param n Total number of genes
#' @param g Number of components
#' @param d Number of data types
#' @param modelname The model name, default is IMIX_ind
#' @return AIC/BIC values of the target model
#' @export

model_selection=function(loglik, #Full log likelihood, result output: `Full MaxLogLik final`
                         n, #Total number of genes
                         g=4, #Number of components
                         d=2, #Number of data types
                         modelname=c("IMIX_ind","IMIX_ind_unrestrict","IMIX_cor_twostep","IMIX_cor","IMIX_cor_restrict") #The model name, default is IMIX_ind
                         ){

  modelname <- match.arg(modelname)
  if(modelname=="IMIX_ind_unrestrict"){
    nparam = g*d + g*d + g-1
  } else if(modelname=="IMIX_ind"){
    nparam = 2*d + 2*d + g-1 #this is independent model: g-1 is the proportion pi
  } else if(modelname=="IMIX_cor_twostep") {
    nparam = g * d*(d+1)/2 + (g - 1) #this is cor fixed mu model.
  }
  else if(modelname=="IMIX_cor") {
    nparam = g * d*(d+1)/2 + g*d + (g - 1) # this is unconstrain cor model: g is number of components, d is dimention, (g - 1) is the proportion pi: (g-1); g*d is number of means; g * d*(d+1)/2 is number of cov parameters.
  }
  else {
    nparam = g * d*(d+1)/2 + 2*d + (g - 1)
  }


  bic=-2*loglik + nparam*log(n)
  aic=-2*loglik + nparam*2
  res=c(aic,bic)
  names(res)=c("AIC","BIC")
  return(res)
}
