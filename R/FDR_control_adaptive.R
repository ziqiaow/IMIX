#' @title The Adaptive Procedure for Across-Data-Type FDR Control
#' @description The adaptive procedure for across-data-type FDR control based on the output from IMIX models, this can be directly performed by IMIX function, however, if the user is interested in other mixture models, alpha level or combinations of components, this function would be useful.
#'
#' @param lfdr Local FDR for each gene of the mixture model results for one component or a combination of components
#' @param alpha Prespecified FDR control level
#' @return The estimated mFDR for the target component or component combinaitons and whether the genes is classified in this component/combination after FDR control at alpha level, 1 is yes, 0 is no.
#' \item{significant_genes_with_FDRcontrol}{The output of each gene ordered by the components based on FDR control and within each component ordered by the local FDR, "localFDR" is 1-posterior probability of each gene in the component based on the maximum posterior probability, "class_withoutFDRcontrol" is the classified component based on maximum posterior probability, "class_FDRcontrol" is the classified component based on the across-data-type FDR control at alpha level}
#' \item{estimatedFDR}{The estimated marginal FDR value for each component starting from component 2 (component 1 is the global null)}
#' \item{alpha}{Prespecified nominal level for the across-data-type FDR control}
#' @export
#' @references
#' Wang, Ziqiao, and Peng Wei. 2020. “IMIX: A Multivariate Mixture Model Approach to Integrative Analysis of Multiple Types of Omics Data.” BioRxiv. Cold Spring Harbor Laboratory. \url{https://doi.org/10.1101/2020.06.23.167312}.
#' @examples 
#' # First load the data
#' data("data_p")
#' 
#' # Specify inititial values (this step could be omitted)
#' mu_input <- c(0,3,0,3)
#' sigma_input <- rep(1,4)
#' p_input <- rep(0.5,4)
#' test1 <- IMIX(data_input = data_p,data_type = "p",mu_ini = mu_input,sigma_ini = sigma_input,
#' p_ini = p_input,alpha = 0.1,model_selection_method = "AIC")
#' 
#' # Check the selected model based on AIC value
#' test1$`Selected Model`
#' 
#' # Below is an example for data example 1 in controlling the FDR at 0.2 for component 2 & component 4. 
#' # First calculate the local FDR for component 2 & component 4:
#' lfdr_ge_combined <- 1 - (test1$IMIX_cor_twostep$`posterior prob`[,2] + 
#' test1$IMIX_cor_twostep$`posterior prob`[,4])  # Class 2: (ge+,cnv-); class 4: (ge+,cnv+)
#' names(lfdr_ge_combined) <- rownames(test1$IMIX_cor_twostep$`posterior prob`)
#' 
#' # Perform the across-data-type FDR control for component 2 & component 4 at alpha level 0.2
#' fdr_control1 <- FDR_control_adaptive(lfdr = lfdr_ge_combined, alpha = 0.2)
#' 



FDR_control_adaptive=function(lfdr, #Local FDR for each gene of the mixture model results for one component or a combination of components
                              alpha #Prespecified FDR control level
                              ){

  m = length(lfdr)
  if (is.null(names(lfdr)) == TRUE) {
    names(lfdr) = paste0("gene", 1:m)
  }
  
  lfdr_ordered = lfdr[order(lfdr)]
  sum_lfdr = 0
  i = 1
  while (sum_lfdr[i] <= alpha) {
    i = i + 1
    sum_lfdr[i] = mean(lfdr_ordered[1:(i - 1)])
  }
  k = i - 2
  if (k == 0) {
    mFDR = 0
    name_rej = NA
    
  } else {
    mFDR = sum_lfdr[k + 1]
    name_rej = names(lfdr_ordered)[1:k]
  }
  
  pred_group = rep(0, m)
  pred_group[match(name_rej, names(lfdr))] = 1
  names(pred_group) = names(lfdr)
  res = list("significant_genes_with_FDRcontrol" = pred_group, "estimatedFDR" = mFDR, "alpha" = alpha)
  return(res)
  
}





#' @title The Adaptive Procedure for Across-Data-Type FDR Control for IMIX Output
#' @description The adaptive procedure for across-data-type FDR control based on the output from IMIX models, this can be directly performed by IMIX function, however, if the user is interested in other alpha levels, this function would be useful to avoid rerun the IMIX().
#'
#' @param imix_output The result output from IMIX() function, result controlled at alpha level only for one component each time.
#' @param model The target model among "IMIX_ind","IMIX_cor_twostep","IMIX_cor_restrict", and "IMIX_cor". Default is IMIX_ind.
#' @param alpha Prespecified FDR control level.
#' @return The estimated mFDR for the target component and classify the genes in each component after FDR control at alpha level.
#' \item{significant_genes_with_FDRcontrol}{The output of each gene ordered by the components based on FDR control and within each component ordered by the local FDR, "localFDR" is 1-posterior probability of each gene in the component based on the maximum posterior probability, "class_withoutFDRcontrol" is the classified component based on maximum posterior probability, "class_FDRcontrol" is the classified component based on the across-data-type FDR control at alpha level}
#' \item{estimatedFDR}{The estimated marginal FDR value for each component starting from component 2 (component 1 is the global null)}
#' \item{alpha}{Prespecified nominal level for the across-data-type FDR control}
#' @export
#' @references
#' Wang, Ziqiao, and Peng Wei. 2020. “IMIX: A Multivariate Mixture Model Approach to Integrative Analysis of Multiple Types of Omics Data.” BioRxiv. Cold Spring Harbor Laboratory. \url{https://doi.org/10.1101/2020.06.23.167312}.
#' @examples 
#' # First generate the data
#' library(MASS)
#' N <- 1000
#' truelabel <- sample(1:8,prob = rep(0.125, 8),size = N,replace = TRUE)
#' mu1 <- c(0, 5);mu2 <- c(0, 5);mu3 <- c(0, 5)
#' mu1_mv <- c(mu1[1], mu2[1], mu3[1]);mu2_mv <- c(mu1[2], mu2[1], mu3[1]);
#' mu3_mv <- c(mu1[1], mu2[2], mu3[1]);mu4_mv <- c(mu1[1], mu2[1], mu3[2]);
#' mu5_mv <- c(mu1[2], mu2[2], mu3[1]);mu6_mv <- c(mu1[2], mu2[1], mu3[2])
#' mu7_mv <- c(mu1[1], mu2[2], mu3[2]);mu8_mv <- c(mu1[2], mu2[2], mu3[2])
#' cov_sim <- list()
#' for (i in 1:8) {
#'   cov_sim[[i]] <- diag(3)
#'   }
#' data_z <- array(0, c(N, 3))
#' data_z[which(truelabel == 1),] <- mvrnorm(n = length(which(truelabel == 1)),
#' mu = mu1_mv,Sigma = cov_sim[[1]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 2),] <- mvrnorm(n = length(which(truelabel == 2)),
#' mu = mu2_mv,Sigma = cov_sim[[2]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 3),] <- mvrnorm(n = length(which(truelabel == 3)),
#' mu = mu3_mv,Sigma = cov_sim[[3]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 4),] <- mvrnorm(n = length(which(truelabel == 4)),
#' mu = mu4_mv,Sigma = cov_sim[[4]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 5),] <- mvrnorm(n = length(which(truelabel == 5)),
#' mu = mu5_mv,Sigma = cov_sim[[5]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 6),] <- mvrnorm(n = length(which(truelabel == 6)),
#' mu = mu6_mv,Sigma = cov_sim[[6]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 7),] <- mvrnorm(n = length(which(truelabel == 7)),
#' mu = mu7_mv,Sigma = cov_sim[[7]],tol = 1e-6,empirical = FALSE)
#' data_z[which(truelabel == 8),] <- mvrnorm(n = length(which(truelabel == 8)),
#' mu = mu8_mv,Sigma = cov_sim[[8]],tol = 1e-6,empirical = FALSE)
#' rownames(data_z) <- paste0("gene", 1:N)
#' colnames(data_z) <- c("z.methy", "z.ge", "z.cnv")
#' dim(data_z)
#' 
#' # Fit the model
#' test2 <- IMIX(data_input = data_z,data_type = "z",alpha = 0.05,verbose = TRUE)
#' 
#' # Adaptive FDR control at alpha 0.2 for IMIX_cor model
#' fdr_control2 <- FDR_control_adaptive_imix(imix_output = test2, model = "IMIX_cor", 
#' alpha = 0.2) 
#' 


FDR_control_adaptive_imix=function(imix_output, #The result output from IMIX() function, result controlled at alpha level only for one component each time
                                   model=c("IMIX_ind","IMIX_cor_twostep","IMIX_cor_restrict","IMIX_cor"), #The target model, default is IMIX_ind
                                   alpha #Prespecified FDR control level
                                   ){
  
  #name = imix_output$`Selected Model`
  name <- match.arg(model)
  best_model = unlist(imix_output[name], recursive = FALSE, use.names = TRUE)
  names(best_model) = gsub(paste0(name, "."), "", names(best_model))
  g = dim(best_model$`posterior prob`)[2]
  
  #########################################
  #Classes based on posterior probability and the corresponding local FDR
  #########################################
  class_before_controlFDR = apply(best_model$`posterior prob`, 1, which.max)
  localFDR_allgenes = 1 - apply(best_model$`posterior prob`, 1, max)
  
  
  #########################################
  #Adaptive procedure for across-data-type FDR control
  #########################################
  pred_group_adaptive = list()
  pred_group_adaptive_mFDR = list()
  pred_group_adaptive_twoclass = list()
  for (comp in 1:(g - 1)) {
    pred_fdr_twoclass = 1 - best_model$`posterior prob`[, (comp + 1)]
    names(pred_fdr_twoclass) = rownames(best_model$`posterior prob`)
    pred_group_adaptive[[comp]] = FDR_control_adaptive(lfdr = pred_fdr_twoclass, alpha =
                                                         alpha)
    pred_group_adaptive_mFDR[[comp]] = pred_group_adaptive[[comp]][[2]]
    pred_group_adaptive_twoclass[[comp]] = pred_group_adaptive[[comp]][[1]]
  }
  names(pred_group_adaptive_mFDR) = paste0("estimated_mFDR_comp", 2:g)
  
  sig_genes_all = array(0, c(dim(best_model$`posterior prob`)[1], 3))
  
  sig_genes_all[, 1] = localFDR_allgenes
  sig_genes_all[, 2] = class_before_controlFDR
  sig_genes_all[, 3] = 1
  
  rownames(sig_genes_all) = rownames(best_model$`posterior prob`)
  colnames(sig_genes_all) = c("localFDR", "class_withoutFDRcontrol", "class_FDRcontrol")
  
  for (comp in 1:(g - 1)) {
    sig_genes_all[which(pred_group_adaptive_twoclass[[comp]] == 1  & sig_genes_all[,2] == (comp+1) ), 3] = comp + 1
  }
  
  sig_genes_all = data.frame(sig_genes_all)
  sig_genes_all = sig_genes_all[order(-sig_genes_all$class_FDRcontrol, sig_genes_all$localFDR), ]
  
  res = list(
    "significant_genes_with_FDRcontrol" = sig_genes_all,
    "estimatedFDR" = pred_group_adaptive_mFDR,
    "alpha" = alpha
  )
  
  return(res)
  
}
