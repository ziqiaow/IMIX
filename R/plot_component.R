#' @title Plot the AIC or BIC Values for Model Selection
#' @description Plot the result output of model selection for components based on AIC and BIC values in IMIX
#'
#' @import ggplot2
#' @param res_select Result output from function model_selection_component()
#' @param type Which information criteria to use for plot
#' @return Plot for the model selection of components
#' @export
#' @references
#' Wang, Ziqiao, and Peng Wei. 2020. “IMIX: A Multivariate Mixture Model Approach to Integrative Analysis of Multiple Types of Omics Data.” BioRxiv. Cold Spring Harbor Laboratory. \url{https://doi.org/10.1101/2020.06.23.167312}.
#' @examples
#' \dontrun{
#' # First load the data
#' data("data_p")
#' 
#' # Perform model selections on the data
#' select_comp1 <- model_selection_component(data_p, data_type = "p", seed = 20)
#' 
#' # Make a plot for BIC values
#' plot_component(select_comp1, type = "BIC")
#' }


plot_component=function(res_select, # Result output from function model_selection_component()
                        type=c("AIC","BIC") # Which information criteria to use for plot
                        ){
  type <- match.arg(type)
  
  df=rbind(data.frame(res_select$`AIC/BIC`$IMIX_ind_unrestrict), data.frame(res_select$`AIC/BIC`$IMIX_cor_twostep),data.frame(res_select$`AIC/BIC`$IMIX_cor))
  g=dim(res_select$`AIC/BIC`$IMIX_ind_unrestrict)[1]
  df$component=rep(c(1:g),3)
  df$Model=rep(c("IMIX_ind_unrestrict","IMIX_cor_twostep","IMIX_cor"),each=g)
  df$AIC=as.numeric(as.character(df$AIC))
  df$BIC=as.numeric(as.character(df$BIC))
  
  if(type=="AIC"){
  p1<-ggplot(df, aes_string(x="component", y="AIC", group="Model")) +
    geom_line(aes_string(color="Model"))+
    geom_point(aes_string(color="Model")) + theme(legend.position="bottom",plot.title = element_text(size=10))+
    geom_point(data=df[which.min(df$AIC),], aes_string(x="component", y="AIC"), size=2,shape=17)
  
  } else {
    p1<-ggplot(df, aes_string(x="component", y="BIC", group="Model")) +
      geom_line(aes_string(color="Model"))+
      geom_point(aes_string(color="Model")) + theme(legend.position="bottom",plot.title = element_text(size=10))+
      geom_point(data=df[which.min(df$BIC),], aes_string(x="component", y="BIC"), size=2,shape=17)
    
    
  }
  p1
  
}
