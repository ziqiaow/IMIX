#' @title IMIX-Cor-Restrict
#' @description Fitting a correlated multivariate mixture model with restrictions on the mean. Input of summary statistics z scores or p values of two or three data types.
#'
#' @param data_input An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
#' @param data_type Whether the input data is the p values or z scores, default is p value
#' @param mu Initial value for the mean of the independent mixture model distribution. A vector of length 2*d, d is number of data types. Needs to be in a special format that corresponds to the initial value of mu, for example, if d=3, needs to be in the format of (null_1,alternative_1,null_2,alternative_2,null_3,alternative_3).
#' @param cov A list of initial values for the covariance matrices. If there are three data types and 8 components, then the initial is a list of 8 covariance matrices, each matix is 3*3.
#' @param p Initial value for the proportion of the distribution in the Gaussian mixture model
#' @param tol The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed set.seed, default is 10
#' @param verbose Whether to print the full log-likelihood for each iteration, default is FALSE
#' @return A list of the results of IMIX-cor-restrict
#' \item{posterior prob}{Posterior probability of each gene for each component}
#' \item{Full LogLik all}{Full log-likelihood of each iteration}
#' \item{Full MaxLogLik final}{The final log-likelihood of the converged model}
#' \item{iterations}{Number of iterations run}
#' \item{pi}{Estimated proportion of each component, sum to 1}
#' \item{mu}{Estimated mean for the null and alternative of each data type: for two data types (mu10,mu11,mu20,mu21), three data types (mu10,mu11,mu20,mu21,mu30,mu31), mui0 is the null for data type i, mui1 is the alternative for data type i.}
#' \item{cov}{A list of estimated variance-covariance matrix of each component}
#'   
#' @importFrom stats qnorm
#' @importFrom utils tail
#' @export
#' @references
#' Wang, Ziqiao, and Peng Wei. 2020. “IMIX: A Multivariate Mixture Model Approach to Integrative Analysis of Multiple Types of Omics Data.” BioRxiv. Cold Spring Harbor Laboratory. \url{https://doi.org/10.1101/2020.06.23.167312}.


#Only specifies the mu, let the sigmas be unconstrained
IMIX_cor_restrict=function(data_input, #An n x d data frame or matrix of the summary statistics z score or p value, n is the nubmer of genes, d is the number of data types. Each row is a gene, each column is a data type.
                           data_type=c("p","z"), #Whether the input data is the p values or z scores, default is p value
                           mu, #Initial value for the mean of the independent mixture model distribution. A vector of length 2*d, d is number of data types. Needs to be in a special format that corresponds to the initial value of mu, for example, if d=3, needs to be in the format of (null_1,alternative_1,null_2,alternative_2,null_3,alternative_3).
                           cov, #A list of initial values for the covariance matrices. If there are three data types and 8 components, then the initial is a list of 8 covariance matrices, each matix is 3*3.
                           p, #Initial value for the proportion of the distribution in the Gaussian mixture model. A vector of length 2^d, d is the number of data types.
                           tol=1e-6, #The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
                           maxiter=1000, #The maximum number of iteration, default is 1000
                           seed=10,#set.seed, default is 10
                           verbose=FALSE #Whether to print the full log-likelihood for each iteration, default is FALSE
){

  data_type <- match.arg(data_type)
  if(data_type=="p"){data_input=apply(data_input,2,function(x) stats::qnorm(x,lower.tail=F))}
  

  n_data=dim(data_input)[2]
  if(length(cov)!=2^n_data | length(mu)!=2*n_data | length(p)!=(2^n_data) | dim(cov[[1]])[1]!=n_data | dim(cov[[1]])[2]!=n_data ) {cat(crayon::red("Error: The dimensions of initial values don't match with each other!")); return(1)}
  
  
  for(i in 2:(2^n_data)) {
    if(dim(cov[[i]])[1]!=dim(cov[[i]])[2] | dim(cov[[i]])[1]!=n_data )  {cat(crayon::red("Error: The dimensions of initial values don't match with each other!")); return(1)}
  }
  
  set.seed(seed)
  
  cat(crayon::cyan$bold("Start IMIX-cor-restrict procedure!\n"))

  # modified sum only considers finite values
  sum.finite <- function(x) {
    sum(x[is.finite(x)])
  }


  if(dim(data_input)[2]==2){
  x=data_input[,1]
  y=data_input[,2]

  mu10=mu[1]
  mu11=mu[2]
  mu20=mu[3]
  mu21=mu[4]

  cov1=cov[[1]]
  cov2=cov[[2]]
  cov3=cov[[3]]
  cov4=cov[[4]]

  pi1=p[1]
  pi2=p[2]
  pi3=p[3]
  pi4=p[4]

  Q <- 0
  # starting value of expected value of the log likelihood
  Q[2] <- sum.finite(log(pi1)+log(mvtnorm::dmvnorm(data_input, c(mu10,mu20), cov1)))+
    sum.finite(log(pi2)+log(mvtnorm::dmvnorm(data_input, c(mu11,mu20), cov2)))+
    sum.finite(log(pi3)+log(mvtnorm::dmvnorm(data_input, c(mu10,mu21), cov3)))+
    sum.finite(log(pi4)+log(mvtnorm::dmvnorm(data_input, c(mu11,mu21), cov4)))

  if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",1,": loglik=",Q[2],"\n")))

  k <- 2

  while (abs(Q[k]-Q[k-1])>=tol & k<=maxiter) {

    # E step
    comp1 <- pi1 * mvtnorm::dmvnorm(data_input, c(mu10,mu20), cov1)
    comp2 <- pi2 * mvtnorm::dmvnorm(data_input, c(mu11,mu20), cov2)
    comp3 <- pi3 * mvtnorm::dmvnorm(data_input, c(mu10,mu21), cov3)
    comp4 <- pi4 * mvtnorm::dmvnorm(data_input, c(mu11,mu21), cov4)


    comp.sum <- comp1 + comp2+comp3+comp4

    p1 <- comp1/comp.sum
    p2 <- comp2/comp.sum
    p3 <- comp3/comp.sum
    p4 <- comp4/comp.sum

    # M step
    pi1 <- sum(p1) / length(x)
    pi2 <- sum(p2) / length(x)
    pi3 <- sum(p3) / length(x)
    pi4 <- sum(p4) / length(x)


    a1=solve(cov1)
    a2=solve(cov2)
    a3=solve(cov3)
    a4=solve(cov4)


    mu10 <- sum( p1*( a1[1,1]*x+a1[1,2]*(y-mu20) ) + p3*( a3[1,1]*x+a3[1,2]*(y-mu21) ) ) / sum(p1*a1[1,1]+p3*a3[1,1]);
    mu20 <- sum( p1*( a1[2,2]*y+a1[1,2]*(x-mu10) ) + p2*( a2[2,2]*y+a2[1,2]*(x-mu11) ) ) / sum(p1*a1[2,2]+p2*a2[2,2]);
    mu11 <- sum( p2*( a2[1,1]*x+a2[1,2]*(y-mu20) ) + p4*( a4[1,1]*x+a4[1,2]*(y-mu21) ) ) / sum(p2*a2[1,1]+p4*a4[1,1]);
    mu21 <- sum( p3*( a3[2,2]*y+a3[1,2]*(x-mu10) ) + p4*( a4[2,2]*y+a4[1,2]*(x-mu11) ) ) / sum(p3*a3[2,2]+p4*a4[2,2]);

    cov1 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu10,mu20)) %*%  as.matrix(t(as.vector(a) - c(mu10,mu20)))) %*% p1 ,nrow=2) / sum(p1);
    cov2 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu11,mu20)) %*%  as.matrix(t(as.vector(a) - c(mu11,mu20)))) %*% p2,nrow=2) / sum(p2);
    cov3 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu10,mu21)) %*%  as.matrix(t(as.vector(a) - c(mu10,mu21)))) %*% p3,nrow=2) / sum(p3);
    cov4 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu11,mu21)) %*%  as.matrix(t(as.vector(a) - c(mu11,mu21)))) %*% p4,nrow=2) / sum(p4);

    k <- k + 1
    Q[k] <- sum(log(comp.sum))
    if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",k-1,": loglik=",Q[k],"\n")))

  }

  pred.values=data.frame(p1,p2,p3,p4)
  cov.final=list(cov1,cov2,cov3,cov4)
  mu.final=c(mu10,mu11,mu20,mu21)
  pi.final=c(pi1,pi2,pi3,pi4)

  } else if (dim(data_input)[2]==3){

    x=data_input[,1]
    y=data_input[,2]
    z=data_input[,3]
    mu10=mu[1]
    mu11=mu[2]
    mu20=mu[3]
    mu21=mu[4]
    mu30=mu[5]
    mu31=mu[6]
    cov1=cov[[1]]
    cov2=cov[[2]]
    cov3=cov[[3]]
    cov4=cov[[4]]
    cov5=cov[[5]]
    cov6=cov[[6]]
    cov7=cov[[7]]
    cov8=cov[[8]]
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
    Q[2] <- sum.finite(log(pi1)+log(mvtnorm::dmvnorm(data_input, c(mu10,mu20,mu30), cov1)))+
      sum.finite(log(pi2)+log(mvtnorm::dmvnorm(data_input, c(mu11,mu20,mu30), cov2)))+
      sum.finite(log(pi3)+log(mvtnorm::dmvnorm(data_input, c(mu10,mu21,mu30), cov3)))+
      sum.finite(log(pi4)+log(mvtnorm::dmvnorm(data_input, c(mu10,mu20,mu31), cov4)))+
      sum.finite(log(pi5)+log(mvtnorm::dmvnorm(data_input, c(mu11,mu21,mu30), cov5)))+
      sum.finite(log(pi6)+log(mvtnorm::dmvnorm(data_input, c(mu11,mu20,mu31), cov6)))+
      sum.finite(log(pi7)+log(mvtnorm::dmvnorm(data_input, c(mu10,mu21,mu31), cov7)))+
      sum.finite(log(pi8)+log(mvtnorm::dmvnorm(data_input, c(mu11,mu21,mu31), cov8)))
    if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",1,": loglik=",Q[2],"\n")))

    k <- 2

    while (abs(Q[k]-Q[k-1])>=tol & k<=maxiter) {

      # E step
      comp1 <- pi1 * mvtnorm::dmvnorm(data_input, c(mu10,mu20,mu30), cov1)
      comp2 <- pi2 * mvtnorm::dmvnorm(data_input, c(mu11,mu20,mu30), cov2)
      comp3 <- pi3 * mvtnorm::dmvnorm(data_input, c(mu10,mu21,mu30), cov3)
      comp4 <- pi4 * mvtnorm::dmvnorm(data_input, c(mu10,mu20,mu31), cov4)
      comp5 <- pi5 * mvtnorm::dmvnorm(data_input, c(mu11,mu21,mu30), cov5)
      comp6 <- pi6 * mvtnorm::dmvnorm(data_input, c(mu11,mu20,mu31), cov6)
      comp7 <- pi7 * mvtnorm::dmvnorm(data_input, c(mu10,mu21,mu31), cov7)
      comp8 <- pi8 * mvtnorm::dmvnorm(data_input, c(mu11,mu21,mu31), cov8)

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
      pi1 <- sum(p1) / length(x)
      pi2 <- sum(p2) / length(x)
      pi3 <- sum(p3) / length(x)
      pi4 <- sum(p4) / length(x)
      pi5 <- sum(p5) / length(x)
      pi6 <- sum(p6) / length(x)
      pi7 <- sum(p7) / length(x)
      pi8 <- sum(p8) / length(x)



      a1=solve(cov1)
      a2=solve(cov2)
      a3=solve(cov3)
      a4=solve(cov4)
      a5=solve(cov5)
      a6=solve(cov6)
      a7=solve(cov7)
      a8=solve(cov8)

      mu10 <- sum( p1*( a1[1,1]*x+a1[1,2]*(y-mu20)+a1[1,3]*(z-mu30) ) + p3*( a3[1,1]*x+a3[1,2]*(y-mu21)+a3[1,3]*(z-mu30) ) +
                     p4*( a4[1,1]*x+a4[1,2]*(y-mu20)+a4[1,3]*(z-mu31) ) + p7*( a7[1,1]*x+a7[1,2]*(y-mu21)+a7[1,3]*(z-mu31) ) ) / sum(p1*a1[1,1]+p3*a3[1,1]+p4*a4[1,1]+p7*a7[1,1]);#print(paste0("mu10=",mu10))
      mu20 <- sum( p1*( a1[2,2]*y+a1[1,2]*(x-mu10)+a1[2,3]*(z-mu30) ) + p2*( a2[2,2]*y+a2[1,2]*(x-mu11)+a2[2,3]*(z-mu30) ) +
                     p4*( a4[2,2]*y+a4[1,2]*(x-mu10)+a4[2,3]*(z-mu31) ) + p6*( a6[2,2]*y+a6[1,2]*(x-mu11)+a6[2,3]*(z-mu31) ) ) / sum(p1*a1[2,2]+p2*a2[2,2]+p4*a4[2,2]+p6*a6[2,2]);#print(paste0("mu20=",mu20))
      mu30 <- sum( p1*( a1[3,3]*z+a1[1,3]*(x-mu10)+a1[2,3]*(y-mu20) ) + p2*( a2[3,3]*z+a2[1,3]*(x-mu11)+a2[2,3]*(y-mu20) ) +
                     p3*( a3[3,3]*z+a3[1,3]*(x-mu10)+a3[2,3]*(y-mu21) ) + p5*( a5[3,3]*z+a5[1,3]*(x-mu11)+a5[2,3]*(y-mu21) ) ) / sum(p1*a1[3,3]+p2*a2[3,3]+p3*a3[3,3]+p5*a5[3,3]);#print(paste0("mu30=",mu30))
      mu11 <- sum( p2*( a2[1,1]*x+a2[1,2]*(y-mu20)+a2[1,3]*(z-mu30) ) + p5*( a5[1,1]*x+a5[1,2]*(y-mu21)+a5[1,3]*(z-mu30) ) +
                     p6*( a6[1,1]*x+a6[1,2]*(y-mu20)+a6[1,3]*(z-mu31) ) + p8*( a8[1,1]*x+a8[1,2]*(y-mu21)+a8[1,3]*(z-mu31) ) ) / sum(p2*a2[1,1]+p5*a5[1,1]+p6*a6[1,1]+p8*a8[1,1]);#print(paste0("mu11=",mu11))
      mu21 <- sum( p3*( a3[2,2]*y+a3[1,2]*(x-mu10)+a3[2,3]*(z-mu30) ) + p5*( a5[2,2]*y+a5[1,2]*(x-mu11)+a5[2,3]*(z-mu30) ) +
                     p7*( a7[2,2]*y+a7[1,2]*(x-mu10)+a7[2,3]*(z-mu31) ) + p8*( a8[2,2]*y+a8[1,2]*(x-mu11)+a8[2,3]*(z-mu31) ) ) / sum(p3*a3[2,2]+p5*a5[2,2]+p7*a7[2,2]+p8*a8[2,2]);#print(paste0("mu21=",mu21))
      mu31 <- sum( p4*( a4[3,3]*z+a4[1,3]*(x-mu10)+a4[2,3]*(y-mu20) ) + p6*( a6[3,3]*z+a6[1,3]*(x-mu11)+a6[2,3]*(y-mu20) ) +
                     p7*( a7[3,3]*z+a7[1,3]*(x-mu10)+a7[2,3]*(y-mu21) ) + p8*( a8[3,3]*z+a8[1,3]*(x-mu11)+a8[2,3]*(y-mu21) ) )/ sum(p4*a4[3,3]+p6*a6[3,3]+p7*a7[3,3]+p8*a8[3,3]);#print(paste0("mu31=",mu31))


      cov1 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu10,mu20,mu30)) %*%  as.matrix(t(as.vector(a) - c(mu10,mu20,mu30)))) %*% p1 ,nrow=3) / sum(p1);#print(cov1)
      cov2 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu11,mu20,mu30)) %*%  as.matrix(t(as.vector(a) - c(mu11,mu20,mu30)))) %*% p2,nrow=3) / sum(p2);#print(cov2)
      cov3 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu10,mu21,mu30)) %*%  as.matrix(t(as.vector(a) - c(mu10,mu21,mu30)))) %*% p3,nrow=3) / sum(p3);#print(cov3)
      cov4 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu10,mu20,mu31)) %*%  as.matrix(t(as.vector(a) - c(mu10,mu20,mu31)))) %*% p4,nrow=3) / sum(p4);#print(cov4)
      cov5 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu11,mu21,mu30)) %*%  as.matrix(t(as.vector(a) - c(mu11,mu21,mu30)))) %*% p5,nrow=3) / sum(p5);#print(cov5)
      cov6 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu11,mu20,mu31)) %*%  as.matrix(t(as.vector(a) - c(mu11,mu20,mu31)))) %*% p6,nrow=3) / sum(p6);#print(cov6)
      cov7 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu10,mu21,mu31)) %*%  as.matrix(t(as.vector(a) - c(mu10,mu21,mu31)))) %*% p7,nrow=3) / sum(p7);#print(cov7)
      cov8 <- matrix( apply(data_input,1, function(a) as.matrix(as.vector(a) - c(mu11,mu21,mu31)) %*%  as.matrix(t(as.vector(a) - c(mu11,mu21,mu31)))) %*% p8,nrow=3) / sum(p8);#print(cov8)


      k <- k + 1
      Q[k] <- sum(log(comp.sum))
      if (verbose==TRUE) cat(crayon::yellow(paste0("iter=",k-1,": loglik=",Q[k],"\n")))

    }

    pred.values=data.frame(p1,p2,p3,p4,p5,p6,p7,p8)
    cov.final=list(cov1,cov2,cov3,cov4,cov5,cov6,cov7,cov8)
    pi.final=c(pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8)
    mu.final=c(mu10,mu11,mu20,mu21,mu30,mu31)

  } else  {cat(crayon::red("Error: Function does not support the number of data types!")); return(1)}

  loglik.approx <- utils::tail(Q, n=1)
  colnames(pred.values)=paste0("component",1:(2^n_data))
  rownames(pred.values)=rownames(data_input)
  if (k>maxiter) cat(crayon::red(paste0("Warning: Exceed maximum iteration k=",maxiter,".\n"))) else cat(crayon::cyan$bold("Successfully Done!\n"))
  res <- list('posterior prob'=pred.values,'Full LogLik all'=Q[-1],'Full MaxLogLik final' = loglik.approx, 'iterations' = k-1,'pi'=pi.final,'mu'=mu.final,'cov'=cov.final)
  return(res)
}



