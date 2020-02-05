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





#' @title Component Selection
#' @description Model selection for components based on AIC and BIC information criteria for models in IMIX
#'
#' @import mclust
#' @param data_input An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
#' @param data_type Whether the input data is the p values or z scores, default is p value
#' @param tol The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed set.seed, default is 10
#' @param verbose Whether to print the full log-likelihood for each iteration, default is FALSE
#' @return Selected number of components based on AIC and BIC
#' 
#' @export

model_selection_component=function(data_input, #An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
                                   data_type=c("p","z"), #Whether the input data is the p values or z scores, default is p value
                                   tol=1e-6, #The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
                                   maxiter=1000, #The maximum number of iteration, default is 1000
                                   seed=10, #set.seed, default is 10
                                   verbose=FALSE #Whether to print the full log-likelihood for each iteration, default is FALSE
                                   ){
  
  data_type <- match.arg(data_type)
  if(data_type=="p"){data_input=apply(data_input,2,function(x) qnorm(x,lower.tail=F))}
  
  set.seed(seed)
  cat(crayon::cyan$bold("Start Number of Component Selections!\n"))
  
  test_ind=list()
  test_fixedmu=list()
  test_cor=list() 
  
  
  for(i in 1:(2^dim(data_input)[2])){
    if(i ==1){
      cat(crayon::cyan$bold(paste0("Test for ",i," Component!\n")))
    }else{
      cat(crayon::cyan$bold(paste0("Test for ",i," Components!\n")))
    }
    
    test_ind[[i]]=mclust::Mclust(data_input,modelNames = "VVI",G=i)
    mu_vec_ind=list()
    for (j in 1:i){
      mu_vec_ind[[j]]=as.numeric(test_ind[[i]]$parameters$mean[,j])
    }
    cov=lapply(seq(dim(test_ind[[i]]$parameters$variance$sigma)[3]),function(x) test_ind[[i]]$parameters$variance$sigma[,,x])
    
    test_fixedmu[[i]]=IMIX_cor_twostep(data_input=data_input,data_type="z",g=i,mu_vec=mu_vec_ind,cov=cov,p=test_ind[[i]]$parameters$pro,seed=seed,tol=tol,maxiter=maxiter,verbose = verbose)
    test_cor[[i]]=IMIX_cor(data_input=data_input,data_type="z",g=i,mu_vec=mu_vec_ind,cov=cov,p=test_ind[[i]]$parameters$pro,seed=seed,tol=tol,maxiter=maxiter,verbose = verbose)
  }
  
  
  
  
  res_test_ind=array(0,c((2^dim(data_input)[2]),3))
  res_test_fixedmu=array(0,c((2^dim(data_input)[2]),3))
  res_test_cor=array(0,c((2^dim(data_input)[2]),3))
  res_test_ind[,1]=res_test_fixedmu[,1]=res_test_cor[,1]=paste0("component",1:(2^dim(data_input)[2]))
  colnames(res_test_ind)=colnames(res_test_fixedmu)=colnames(res_test_cor)=c("component","AIC","BIC")
  for(i in 1:(2^dim(data_input)[2])){
    res_test_ind[i,2:3]=model_selection(test_ind[[i]]$loglik,dim(data_input)[1],g=i,modelname="IMIX_ind_unrestrict")
    res_test_fixedmu[i,2:3]=model_selection(test_fixedmu[[i]]$`Full MaxLogLik final`,dim(data_input)[1],g=test_fixedmu[[i]]$g,modelname="IMIX_cor_twostep")
    res_test_cor[i,2:3]=model_selection(test_cor[[i]]$`Full MaxLogLik final`,dim(data_input)[1],g=test_cor[[i]]$g,modelname="IMIX_cor")
    
  }
  
  
  res_list=list(res_test_ind,res_test_fixedmu,res_test_cor)
  names(res_list)=c("IMIX_ind_unrestrict","IMIX_cor_twostep","IMIX_cor")
  df = as.data.frame(do.call(rbind, lapply(res_list, unlist)))
  df$AIC=as.numeric(as.character(df$AIC));df$BIC=as.numeric(as.character(df$BIC))
  best_component_aic=df$component[df$AIC==min(df$AIC)]  
  best_component_bic=df$component[df$BIC==min(df$BIC)]  
  
  cat(crayon::cyan$bold("Done!\n"))
  
  res=list("Component_Selected_AIC"=best_component_aic,"Component_Selected_BIC"=best_component_bic,"AIC/BIC"=res_list,"IMIX_ind_unrestrict"=test_ind,"IMIX_cor_twostep"=test_fixedmu,"IMIX_cor"=test_cor)
  
  return(res)
}




