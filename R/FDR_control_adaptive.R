#' @title The adaptive procedure for global FDR control
#' @description The adaptive procedure for global FDR control based on the output from IMIX models, this can be directly performed by IMIX function, however, if you are interested in other mixture models, alpha level or combinations of components, this function would be useful.
#'
#' @param lfdr Local FDR for each gene of the mixture model results for one component or a combination of components
#' @param alpha Prespecified FDR control level
#' @return The estimated mFDR for the target component or component combinaitons and whether the genes is classified in this component/combination after FDR control at alpha level, 1 is yes, 0 is no.
#' @export


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





#' @title The adaptive procedure for global FDR control for IMIX output
#' @description The adaptive procedure for global FDR control based on the output from IMIX models, this can be directly performed by IMIX function, however, if you are interested in other alpha levels, this function would be useful to avoid rerun the IMIX().
#'
#' @param imix_output The result output from IMIX() function, result controled at alpha level only for one component each time
#' @param alpha Prespecified FDR control level
#' @return The estimated mFDR for the target component and classify the genes in each component after FDR control at alpha level.
#' @export

FDR_control_adaptive_imix=function(imix_output, #The result output from IMIX() function, result controled at alpha level only for one component each time
                                   alpha #Prespecified FDR control level
                                   ){
  
  name = imix_output$`Selected Model`
  best_model = unlist(imix_output[name], recursive = FALSE, use.names = TRUE)
  names(best_model) = gsub(paste0(name, "."), "", names(best_model))
  g = dim(best_model$`posterior prob`)[2]
  
  #########################################
  #Classes based on posterior probability and the corresponding local FDR
  #########################################
  class_before_controlFDR = apply(best_model$`posterior prob`, 1, which.max)
  localFDR_allgenes = 1 - apply(best_model$`posterior prob`, 1, max)
  
  
  #########################################
  #Adaptive procedure for global FDR control
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
    sig_genes_all[which(pred_group_adaptive_twoclass[[comp]] == 1), 3] = comp + 1
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
