#' @title IMIX-Cor-Twostep
#' @description Fitting a correlated multivariate mixture model with fixed mean from estimated parameters of IMIX-ind. Input of summary statistics z scores or p values of two or three data types.
#'
#' @param data_input An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
#' @param data_type Whether the input data is the p values or z scores, default is p value
#' @param g The number of components, default is 8 for three data types
#' @param mu_vec Input of the mean value output from IMIX-Ind result, a list of the mean vectors for each component.
#' @param cov A list of initial values for the covariance matrices. If there are three data types and 8 components, then the initial is a list of 8 covariance matrices, each matix is 3*3.
#' @param p Initial value for the proportion of the distribution in the Gaussian mixture model
#' @param tol The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed set.seed, default is 10
#' @param verbose Whether to print the full log-likelihood for each iteration, default is FALSE
#' @return A list of the results of IMIX-cor-twostep
#' \item{posterior prob}{Posterior probability matrix of each gene for each component}
#' \item{Full LogLik all}{Full log-likelihood of each iteration}
#' \item{Full MaxLogLik final}{The final log-likelihood of the converged model}
#' \item{iterations}{Number of iterations run}
#' \item{pi}{Estimated proportion of each component, sum to 1}
#' \item{mu}{A list of mean vectors of each component for each data type, this is the prespecified mean}
#' \item{cov}{A list of estimated variance-covariance matrix of each component}
#' \item{g}{Number of components}
#' 
#' @importFrom stats qnorm
#' @importFrom utils tail
#' @export
#' @references
#' Ziqiao Wang and Peng Wei. 2020. “IMIX: a multivariate mixture model approach to association analysis through multi-omics data integration.” Bioinformatics. \url{https://doi.org/10.1093/bioinformatics/btaa1001}.

IMIX_cor_twostep=function(data_input, #An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
                          data_type=c("p","z"), #Whether the input data is the p values or z scores, default is p value
                          g=8, #The number of components, default is 8 for three data types
                          mu_vec, #Input of the mean value output from IMIX-Ind result, a list of the mean vectors for each component.
                          cov, #A list of initial values for the covariance matrices, if there are three data types and 8 components, then the initial is a list of 8 covariance matrices, each is 3*3.
                          p, #Initial value for the proportion of the distribution in the Gaussian mixture model
                          tol=1e-6, #The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
                          maxiter=1000, #The maximum number of iteration, default is 1000
                          seed=10,#set.seed, default is 10
                          verbose=FALSE #Whether to print the full log-likelihood for each iteration, default is FALSE
){
 
  data_type <- match.arg(data_type)
  if(data_type=="p"){data_input=apply(data_input,2,function(x) stats::qnorm(x,lower.tail=F))}
  
  n_data=dim(data_input)[2]
  if(length(cov)!=g | length(mu_vec)!=g | length(mu_vec[[1]])!=n_data | dim(cov[[1]])[1]!=n_data | dim(cov[[1]])[2]!=n_data | length(p)!=g | g>(2^n_data) ) {cat(crayon::red("Error: The dimensions of initial values don't match with each other!")); return(1)}
  if(g>1){
  for(i in 2:g) {
    if(dim(cov[[i]])[1]!=dim(cov[[i]])[2] | dim(cov[[i]])[1]!=n_data | length(mu_vec[[i]])!=length(mu_vec[[1]]))  {cat(crayon::red("Error: The dimensions of initial values don't match with each other!")); return(1)}
  }
  }
  set.seed(seed)
  

  cat(crayon::cyan$bold("Start IMIX-cor-twostep procedure!\n"))
  # modified sum only considers finite values
  # sum.finite <- function(x) {
  #    sum(x[is.finite(x)])
  #  }


  N=dim(data_input)[1]
  Q <- 0

  i=1
  if(g==1){
    Q[2] <- sum(log(p[i])+log(mvtnorm::dmvnorm(data_input, mu_vec[[i]], cov[[i]])))
  } else {
    Q[2] <- sum(log(p[i])+log(mvtnorm::dmvnorm(data_input, mu_vec[[i]], cov[[i]])))
    for(i in 2:g){
      Q[2] <- Q[2] + sum(log(p[i])+log(mvtnorm::dmvnorm(data_input, mu_vec[[i]], cov[[i]])))

    }
  }
  if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",1,": loglik=",Q[2],"\n")))

  k <- 2

  while (abs(Q[k]-Q[k-1])>=tol & k<=maxiter) {

    comp=array(0,c(g,N))
    for(i in 1:g){
      comp[i,] <- p[i] * mvtnorm::dmvnorm(data_input, mu_vec[[i]], cov[[i]])
    }

    comp.sum <- apply(comp,2,sum)

    proportion <- apply(comp,1, function(x) x / comp.sum)

    # M step
    p=0

    for(i in 1:g){
      p[i] <- sum(proportion[,i]) / dim(data_input)[1]
    }


    cov=list()
    for(i in 1:g){
      cov[[i]] <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - mu_vec[[i]]) %*%  as.matrix(t(as.vector(a) - mu_vec[[i]]))) %*% proportion[,i] ,nrow=dim(data_input)[2]) / sum(proportion[,i])
    }

    k <- k + 1
    Q[k] <- sum(log(comp.sum))
    if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",k-1,": loglik=",Q[k],"\n")))

  }

  loglik.approx <- utils::tail(Q, n=1)
  pred.values=proportion
  colnames(pred.values)=paste0("component",1:g)
  rownames(pred.values)=rownames(data_input)
  if (k>maxiter) cat(crayon::red(paste0("Warning: Exceed maximum iteration k=",maxiter,".\n"))) else cat(crayon::cyan$bold("Successfully Done!\n"))
  res <- list('posterior prob'=pred.values,'Full LogLik all'=Q[-1],'Full MaxLogLik final' = loglik.approx,'iterations' = k-1,'pi'=p,'mu'=mu_vec,'cov'=cov,'g'=g)
  return(res)

}





