#' @title IMIX-ind
#' @description Fitting an independent mixture model with restrictions on mean and variance. Input of summary statistics z scores or p values of two or three data types.
#'
#' @param data_input An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
#' @param data_type Whether the input data is the p values or z scores, default is p value
#' @param mu Initial value for the mean of each component of the independent mixture model distribution. A vector of length 2*d, d is number of data types. Needs to be in a special format that corresponds to the initial value of mu, for example, if d=3, needs to be in the format of (null_1,alternative_1,null_2,alternative_2,null_3,alternative_3).
#' @param sigma Initial value for the standard deviation of each component of the independent mixture model distribution. A vector of length 2*d, d is number of data types. Needs to be in a special format that corresponds to the initial value of mu, for example, if d=3, needs to be in the format of (null_1,alternative_1,null_2,alternative_2,null_3,alternative_3).
#' @param p Initial value for the proportion of the distribution in the Gaussian mixture model
#' @param tol The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed set.seed, default is 10
#' @param verbose Whether to print the full log-likelihood for each iteration, default is FALSE
#' @return A list of the results of IMIX-ind
#' \item{posterior prob}{Posterior probability matrix of each gene for each component}
#' \item{Full LogLik all}{Full log-likelihood of each iteration}
#' \item{Full MaxLogLik final}{The final log-likelihood of the converged model}
#' \item{iterations}{Number of iterations run}
#' \item{pi}{Estimated proportion of each component, sum to 1}
#' \item{mu}{Estimated mean for the null and alternative of each data type: for two data types (mu10,mu11,mu20,mu21), three data types (mu10,mu11,mu20,mu21,mu30,mu31), mui0 is the null for data type i, mui1 is the alternative for data type i.}
#' \item{sigma}{Estimated standard deviation for the null and alternative of each data type: for two data types (sigma10,sigma11,sigma20,sigma21), three data types (sigma10,sigma11,sigma20,sigma21,sigma30,sigma31), sigmai0 is the null for data type i, sigmai1 is the alternative for data type i.}
#' 
#' @importFrom stats qnorm dnorm
#' @importFrom utils tail
#' @export
#' @references
#' Ziqiao Wang and Peng Wei. 2020. “IMIX: a multivariate mixture model approach to association analysis through multi-omics data integration.” Bioinformatics. \url{https://doi.org/10.1093/bioinformatics/btaa1001}.

IMIX_ind=function(data_input, #An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
                  data_type=c("p","z"), #Whether the input data is the p values or z scores, default is p value
                  mu, #Initial value for the mean of the independent mixture model distribution. A vector of length 2*d, d is number of data types. Needs to be in a special format that corresponds to the initial value of mu, for example, if d=3, needs to be in the format of (null_1,alternative_1,null_2,alternative_2,null_3,alternative_3).
                  sigma, #Initial value for the standard deviation of the independent mixture model distribution. A vector of length 2*d, d is number of data types. Needs to be in a special format that corresponds to the initial value of mu, for example, if d=3, needs to be in the format of (null_1,alternative_1,null_2,alternative_2,null_3,alternative_3).
                  p, #Initial value for the proportion of the distribution in the Gaussian mixture model
                  tol=1e-6, #The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
                  maxiter=1000, #The maximum number of iteration, default is 1000
                  seed=10,#set.seed, default is 10
                  verbose=FALSE #Whether to print the full log-likelihood for each iteration, default is FALSE
                  ){
  
  
  data_type <- match.arg(data_type)
  if(data_type=="p"){data_input=apply(data_input,2,function(x) stats::qnorm(x,lower.tail=F))}
  


  n_data=dim(data_input)[2]
  if(length(sigma)!=2*n_data | length(mu)!=2*n_data | length(p)!=(2^n_data) ) {cat(crayon::red("Error: The dimensions of initial values don't match with each other!")); return(1)}

  
  set.seed(seed)
  
  
  cat(crayon::cyan$bold("Start IMIX-ind procedure!\n"))

  # modified sum only considers finite values
  sum.finite <- function(x) {
    sum(x[is.finite(x)])
  }


  if(dim(data_input)[2]==3){
  x=data_input[,1]
  y=data_input[,2]
  z=data_input[,3]
  mu10=mu[1]
  mu11=mu[2]
  mu20=mu[3]
  mu21=mu[4]
  mu30=mu[5]
  mu31=mu[6]
  sigma10=sigma[1]
  sigma11=sigma[2]
  sigma20=sigma[3]
  sigma21=sigma[4]
  sigma30=sigma[5]
  sigma31=sigma[6]
  pi1=p[1]
  pi2=p[2]
  pi3=p[3]
  pi4=p[4]
  pi5=p[5]
  pi6=p[6]
  pi7=p[7]
  pi8=p[8]

  Q <- 0
  # starting value of expected value of the log likelihood
  Q[2] <- sum.finite(log(pi1)+log(stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu30,sigma30))) +
    sum.finite(log(pi2)+log(stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu30,sigma30)))+
    sum.finite(log(pi3)+log(stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu30,sigma30)))+
    sum.finite(log(pi4)+log(stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu31,sigma31)))+
    sum.finite(log(pi5)+log(stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu30,sigma30)))+
    sum.finite(log(pi6)+log(stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu31,sigma31)))+
    sum.finite(log(pi7)+log(stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu31,sigma31)))+
    sum.finite(log(pi8)+log(stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu31,sigma31)))
  if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",1,": loglik=",Q[2],"\n")))

  k <- 2

  while (abs(Q[k]-Q[k-1])>=tol & k<=maxiter) {
    # E step
    comp1 <- pi1 * stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu30,sigma30)
    comp2 <- pi2 * stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu30,sigma30)
    comp3 <- pi3 * stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu30,sigma30)
    comp4 <- pi4 * stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu31,sigma31)
    comp5 <- pi5 * stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu30,sigma30)
    comp6 <- pi6 * stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu20,sigma20)*stats::dnorm(z,mu31,sigma31)
    comp7 <- pi7 * stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu31,sigma31)
    comp8 <- pi8 * stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu21,sigma21)*stats::dnorm(z,mu31,sigma31)

    comp.sum <- comp1 + comp2+comp3+comp4+comp5+comp6+comp7+comp8

    p1 <- comp1/comp.sum
    p2 <- comp2/comp.sum
    p3 <- comp3/comp.sum
    p4 <- comp4/comp.sum
    p5 <- comp5/comp.sum
    p6 <- comp6/comp.sum
    p7 <- comp7/comp.sum
    p8 <- comp8/comp.sum

    # M step
    pi1 <- sum.finite(p1) / length(x)
    pi2 <- sum.finite(p2) / length(x)
    pi3 <- sum.finite(p3) / length(x)
    pi4 <- sum.finite(p4) / length(x)
    pi5 <- sum.finite(p5) / length(x)
    pi6 <- sum.finite(p6) / length(x)
    pi7 <- sum.finite(p7) / length(x)
    pi8 <- sum.finite(p8) / length(x)

    mu10 <- sum.finite((p1+p3+p4+p7) * x) / sum.finite(p1+p3+p4+p7)
    mu20 <- sum.finite((p1+p2+p4+p6) * y) / sum.finite(p1+p2+p4+p6)
    mu30 <- sum.finite((p1+p2+p3+p5) * z) / sum.finite(p1+p2+p3+p5)
    mu11 <- sum.finite((p2+p5+p6+p8) * x) / sum.finite(p2+p5+p6+p8)
    mu21 <- sum.finite((p3+p5+p7+p8) * y) / sum.finite(p3+p5+p7+p8)
    mu31 <- sum.finite((p4+p6+p7+p8) * z) / sum.finite(p4+p6+p7+p8)


    sigma10 <- sqrt(sum.finite((p1+p3+p4+p7) * (x-mu10)^2) / sum.finite(p1+p3+p4+p7))
    sigma20 <- sqrt(sum.finite((p1+p2+p4+p6) * (y-mu20)^2) / sum.finite(p1+p2+p4+p6))
    sigma30 <- sqrt(sum.finite((p1+p2+p3+p5) * (z-mu30)^2) / sum.finite(p1+p2+p3+p5))
    sigma11 <- sqrt(sum.finite((p2+p5+p6+p8) * (x-mu11)^2) / sum.finite(p2+p5+p6+p8))
    sigma21 <- sqrt(sum.finite((p3+p5+p7+p8) * (y-mu21)^2) / sum.finite(p3+p5+p7+p8))
    sigma31 <- sqrt(sum.finite((p4+p6+p7+p8) * (z-mu31)^2) / sum.finite(p4+p6+p7+p8))


    k <- k + 1
    Q[k] <- sum(log(comp.sum))
    if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",k-1,": loglik=",Q[k],"\n")))
  }
  pred.values=data.frame(p1,p2,p3,p4,p5,p6,p7,p8)
  pi_final=c(pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8);mu_final=c(mu10,mu11,mu20,mu21,mu30,mu31);sigma_final=c(sigma10,sigma11,sigma20,sigma21,sigma30,sigma31)

   } else if (dim(data_input)[2]==2){
    x=data_input[,1]
    y=data_input[,2]
    mu10=mu[1]
    mu11=mu[2]
    mu20=mu[3]
    mu21=mu[4]

    sigma10=sigma[1]
    sigma11=sigma[2]
    sigma20=sigma[3]
    sigma21=sigma[4]

    pi1=p[1]
    pi2=p[2]
    pi3=p[3]
    pi4=p[4]


    Q <- 0

    # starting value of expected value of the log likelihood
    Q[2] <- sum.finite(log(pi1)+log(stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu20,sigma20))) +
      sum.finite(log(pi2)+log(stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu20,sigma20)))+
      sum.finite(log(pi3)+log(stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu21,sigma21)))+
      sum.finite(log(pi4)+log(stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu21,sigma21)))

    if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",1,": loglik=",Q[2],"\n")))

    k <- 2

    while (abs(Q[k]-Q[k-1])>=tol & k<=maxiter) {
      # E step
      comp1 <- pi1 * stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu20,sigma20)
      comp2 <- pi2 * stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu20,sigma20)
      comp3 <- pi3 * stats::dnorm(x, mu10, sigma10)*stats::dnorm(y,mu21,sigma21)
      comp4 <- pi4 * stats::dnorm(x, mu11, sigma11)*stats::dnorm(y,mu21,sigma21)


      comp.sum <- comp1 + comp2+comp3+comp4

      p1 <- comp1/comp.sum
      p2 <- comp2/comp.sum
      p3 <- comp3/comp.sum
      p4 <- comp4/comp.sum


      # M step
      pi1 <- sum.finite(p1) / length(x)
      pi2 <- sum.finite(p2) / length(x)
      pi3 <- sum.finite(p3) / length(x)
      pi4 <- sum.finite(p4) / length(x)


      mu10 <- sum.finite((p1+p3) * x) / sum.finite(p1+p3)
      mu20 <- sum.finite((p1+p2) * y) / sum.finite(p1+p2)

      mu11 <- sum.finite((p2+p4) * x) / sum.finite(p2+p4)
      mu21 <- sum.finite((p3+p4) * y) / sum.finite(p3+p4)



      sigma10 <- sqrt(sum.finite((p1+p3) * (x-mu10)^2) / sum.finite(p1+p3))
      sigma20 <- sqrt(sum.finite((p1+p2) * (y-mu20)^2) / sum.finite(p1+p2))

      sigma11 <- sqrt(sum.finite((p2+p4) * (x-mu11)^2) / sum.finite(p2+p4))
      sigma21 <- sqrt(sum.finite((p3+p4) * (y-mu21)^2) / sum.finite(p3+p4))



      k <- k + 1
      Q[k] <- sum(log(comp.sum))
      if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",k-1,": loglik=",Q[k],"\n")))
    }
    pred.values=data.frame(p1,p2,p3,p4)
    pi_final=c(pi1,pi2,pi3,pi4);mu_final=c(mu10,mu11,mu20,mu21);sigma_final=c(sigma10,sigma11,sigma20,sigma21)
   } else  {cat(crayon::red("Error: Function does not support the number of data types!")); return(1)}

  loglik.approx <- utils::tail(Q, n=1)
  colnames(pred.values)=paste0("component",1:(2^n_data))
  rownames(pred.values)=rownames(data_input)
  if (k>maxiter) cat(crayon::red(paste0("Warning: Exceed maximum iteration k=",maxiter,".\n"))) else cat(crayon::cyan$bold("Successfully Done!\n"))
  res <- list('posterior prob'=pred.values,'Full LogLik all'=Q[-1],'Full MaxLogLik final' = loglik.approx, 'iterations' = k-1,'pi'=pi_final,'mu'=mu_final,'sigma'=sigma_final)
  return(res)
}
