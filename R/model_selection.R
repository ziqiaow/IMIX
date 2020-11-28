#' @title Model Selection
#' @description Model selection for sub-model outputs in IMIX, this step is to calculate the AIC or BIC values for one model
#'
#' @param loglik Full log likelihood, result output from IMIX or a sub model in IMIX: `Full MaxLogLik final`
#' @param n Total number of genes
#' @param g Number of components
#' @param d Number of data types
#' @param modelname The model name, default is IMIX_ind
#' @return AIC/BIC values of the target model
#' 
#' @export
#' @references
#' Wang, Ziqiao, and Peng Wei. 2020. “IMIX: A Multivariate Mixture Model Approach to Integrative Analysis of Multiple Types of Omics Data.” BioRxiv. Cold Spring Harbor Laboratory. \url{https://doi.org/10.1101/2020.06.23.167312}.
#' @examples 
#' # First load the data
#' data("data_p")
#' 
#' # Specify the initial values
#' mu_input <- c(0,3,0,3)
#' sigma_input <- rep(1,4)
#' p_input <- rep(0.5,4)
#' 
#' # Fit the IMIX model
#' test1 <- IMIX(data_input = data_p,data_type = "p",mu_ini = mu_input,sigma_ini = sigma_input,
#' p_ini = p_input,alpha = 0.1,model_selection_method = "AIC")
#' 
#' # Calculate the AIC and BIC values for IMIX_ind with two data types and four components
#' model_selection(test1$IMIX_ind$`Full MaxLogLik final`,
#' n=dim(test1$IMIX_ind$`posterior prob`)[1],g=4,d=2, "IMIX_ind")



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


  BIC=-2*loglik + nparam*log(n)
  AIC=-2*loglik + nparam*2
  res=c(AIC,BIC)
  names(res)=c("AIC","BIC")
  return(res)
}





#' @title Component Selection
#' @description Model selection for components based on AIC and BIC values for models in IMIX
#'
#' @import mclust
#' @param data_input An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
#' @param data_type Whether the input data is the p values or z scores, default is p value
#' @param tol The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed set.seed, default is 10
#' @param verbose Whether to print the full log-likelihood for each iteration, default is FALSE
#' @return Selected number of components based on AIC and BIC
#' \item{Component_Selected_AIC}{Selected number of components by AIC with the smallest AIC value among all components and models}
#' \item{Component_Selected_BIC}{Selected number of components by BIC with the smallest BIC value among all components and models}
#' \item{AIC/BIC}{The AIC and BIC values for all components for IMIX_ind_unrestrict, IMIX_cor_twostep, and IMIX_cor}
#' \item{IMIX_ind_unrestrict}{A list of the IMIX_ind_unrestrict for all components 1,2,...2^d, this step was fitted using R package "Mclust", more details of the output can be found there}
#' \item{IMIX_cor_twostep}{A list of the IMIX_cor_twostep for all components 1,2,...2^d, here, the mean is the estimated value of IMIX_ind_unrestrict}
#' \item{IMIX_cor}{A list of the IMIX_cor_twostep for all components 1,2,...2^d}
#' 
#' @importFrom stats qnorm
#' @export
#' @references
#' Wang, Ziqiao, and Peng Wei. 2020. “IMIX: A Multivariate Mixture Model Approach to Integrative Analysis of Multiple Types of Omics Data.” BioRxiv. Cold Spring Harbor Laboratory. \url{https://doi.org/10.1101/2020.06.23.167312}.
#' 
#' Scrucca, Luca, Michael Fop, T. Brendan Murphy, and Adrian E. Raftery. 2016. “mclust 5: Clustering, Classification and Density Estimation Using Gaussian Finite Mixture Models.” The R Journal 8 (1): 289–317. \url{https://doi.org/10.32614/RJ-2016-021}.
#' @examples
#' 
#' # First load the data
#' data("data_p")
#' 
#' # Perform model selections on the data
#' select_comp1 = model_selection_component(data_p, data_type = "p", seed = 20)


model_selection_component=function(data_input, #An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
                                   data_type=c("p","z"), #Whether the input data is the p values or z scores, default is p value
                                   tol=1e-6, #The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
                                   maxiter=1000, #The maximum number of iteration, default is 1000
                                   seed=10, #set.seed, default is 10
                                   verbose=FALSE #Whether to print the full log-likelihood for each iteration, default is FALSE
                                   ){
  
  data_type <- match.arg(data_type)
  if(data_type=="p"){data_input=apply(data_input,2,function(x) stats::qnorm(x,lower.tail=F))}
  
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
  best_component_AIC=df$component[df$AIC==min(df$AIC)]  
  best_component_BIC=df$component[df$BIC==min(df$BIC)]  
  
  cat(crayon::cyan$bold("Done!\n"))
  
  res=list("Component_Selected_AIC"=best_component_AIC,"Component_Selected_BIC"=best_component_BIC,"AIC/BIC"=res_list,"IMIX_ind_unrestrict"=test_ind,"IMIX_cor_twostep"=test_fixedmu,"IMIX_cor"=test_cor)
  
  return(res)
}




