#' @title The adaptive procedure for global FDR control
#' @description The adaptive procedure for global FDR control based on the output from IMIX models, this can be directly performed by IMIX function, however, if you are interested in other mixture models, this function would be useful.
#'
#' @param lfdr Local FDR for each gene of the mixture model results for one component
#' @param alpha Prespecified FDR control level
#' @return The estimated mFDR for the target component and whether the genes is classified in this component after FDR control at alpha level, 1 is yes, 0 is no.
#' @export
FDR_control_adaptive=function(lfdr, #Local FDR for each gene of the mixture model results for one component
                              alpha #Prespecified FDR control level
                              ){
  m=length(lfdr)
  if(is.null(names(lfdr))==TRUE){
    names(lfdr)=paste0("gene",1:m)}

  lfdr_ordered=lfdr[order(lfdr)]
  sum_lfdr=0
  i=1
  while(sum_lfdr[i]<=alpha){
    i=i+1
    sum_lfdr[i]=mean(lfdr_ordered[1:(i-1)])
  }
  k=i-2
  if(k==0){
    mFDR=0
    name_rej=NA

  } else {
    mFDR=sum_lfdr[k+1]
    name_rej=names(lfdr_ordered)[1:k]}

  pred_group=rep(0,m)
  pred_group[match(name_rej,names(lfdr))]=1
  names(pred_group)=names(lfdr)
  res=list(pred_group,mFDR)
  return(res)
}

