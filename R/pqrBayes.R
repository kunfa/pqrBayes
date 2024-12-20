#' fit a regularized Bayesian quantile varying coefficient model
#' 
#' @keywords models
#' @param g the matrix of predictors (subject to selection) without intercept.
#' @param y the response variable. The current version only supports the continuous response.
#' @param u a vector of effect modifying variable of the quantile varying coefficient model.
#' @param e a matrix of clinical covariates not subject to selection.
#' @param quant the quantile level specified by users. The default value is 0.5.
#' @param iterations the number of MCMC iterations.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B-spline basis.
#' @param robust logical flag. If TRUE, robust methods will be used.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly.
#' @param hyper a named list of hyperparameters.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#'
#' @details The model described in "\code{\link{data}}" is:
#' \deqn{Y_{i}=\sum_{k=1}^{q} E_{ik} \beta_k +\sum_{j=0}^{p}\gamma_j(V_i)X_{ij} +\epsilon_{i},}
#' where \eqn{\beta_k}'s are the regression coefficients for the clinical covariates and \eqn{\gamma_j}'s are the varying coefficients for the intercept and predictors (e.g. genetic factors).
#'
#' 
#' When \{sparse=TRUE\} (default), spike--and--slab priors are adopted. Otherwise, Laplacian shrinkage will be used.
#
#
#' Users can modify the hyper-parameters by providing a named list of hyper-parameters via the argument `hyper'.
#' The list can have the following named components
#' \describe{
#'   \item{a0, b0: }{ shape parameters of the Beta priors (\eqn{\pi^{a_{0}-1}(1-\pi)^{b_{0}-1}}) on \eqn{\pi_{0}}.}
#'   \item{c1, c2: }{ the shape parameter and the rate parameter of the Gamma prior on \eqn{\nu}.}
#' }
#' 
#' Please check the references for more details about the prior distributions.
#' 
#' @return an object of class "pqrBayes" is returned, which is a list with components:
#' \item{posterior}{posterior samples from the MCMC}
#' \item{coefficients}{a list of posterior estimates of coefficients}
#' 
#' @examples
#' data(data)
#' g=data$g
#' y=data$y
#' u=data$u
#' e=data$e
#'
#' ## default method
#' fit1=pqrBayes(g,y,u,e,quant=0.5)
#' fit1
#'
#' \donttest{
#'
#' ## non-sparse
#' sparse=FALSE
#' fit2=pqrBayes(g,y,u,e,quant=0.5,sparse = sparse)
#' fit2
#' }
#' @export
#' @references
#' 
#' Zhou, F., Ren, J., Ma, S. and Wu, C. (2023). The Bayesian regularized quantile varying coefficient model.
#'  {\emph{Computational Statistics & Data Analysis}, 107808} \doi{10.1016/j.csda.2023.107808}
#'  
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, 79(2), 684-694} \doi{10.1111/biom.13670}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \doi{10.1002/sim.8434}


pqrBayes <- function(g, y, u, e=NULL,quant=0.5, iterations=10000, kn=2, degree=2, robust = TRUE,sparse=TRUE, hyper=NULL,debugging=FALSE){
  p = dim(g)[2]
  
  x = cbind(1,g)
  n = length(y); 
  
  
  
  ## basis expansion
  
  d=kn+degree+1
  u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots = as.numeric(stats::quantile(u, u.k))
  pi.u = splines::bs(u, knots=Knots, intercept=TRUE, degree=degree)[,1:(d)]
  
  
  xx = as.data.frame(matrix(0, n, (p+1)*d))
  for(j in 1:(p+1)){
    last = j*d; first = last-d+1
    xx[,first:last] = pi.u*x[,j]
  }
  xx = as.matrix(xx)
  
  if(!is.null(e)){
    q = dim(e)[2]
    xxwe = cbind(xx, e)
    lasso.cv = glmnet::cv.glmnet(xxwe,y,alpha=1,nfolds=5)
    lambda.cv = lasso.cv$lambda.min;
    lasso.fit = glmnet::glmnet(xxwe, y, family="gaussian",alpha=1,nlambda=50)
    coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-1];
    
    hat.m = coeff.array[1:d]      ## coeff for varying intercept
    hat.r = coeff.array[(d+1):((p+1)*d)] ## coeff for varying part
    hat.clin = coeff.array[((p+1)*d+1):dim(xxwe)[2]] ## coeff for clinic covariates
    
    
    xx1=xx[,-(1:d)]
    CLC=cbind(pi.u,e)
    hatAlpha=c(hat.m,hat.clin)
    invSigAlpha0 = diag(10^-3, (d+q))
  }else{
    q = 0
    xxwe = xx
    lasso.cv = glmnet::cv.glmnet(xxwe,y,alpha=1,nfolds=5)
    lambda.cv = lasso.cv$lambda.min;
    lasso.fit = glmnet::glmnet(xxwe, y, family="gaussian",alpha=1,nlambda=50)
    coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-1];
    
    hat.m = coeff.array[1:d]      ## coeff for varying intercept
    hat.r = coeff.array[(d+1):((p+1)*d)] ## coeff for varying part
    
    
    xx1=xx[,-(1:d)]
    CLC=pi.u
    hatAlpha=hat.m
    invSigAlpha0 = diag(10^-3, (d+q))
  }
  
  xi1 = (1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatTau = 1
  hatV = rep(1,n) 
  hatEtaSq = 1
  hatSg = rep(1, p)
  hatPi = 0.5
  invTAUsq.star = rep(0.1, p)
  hat.pi.s = 0.9
  lambda.star = 1
  hat.sigma.sq = 1
  a.star=1
  b.star = 1.5
  alpha = 0.2
  gamma = 0.1
  mu.star=nu.star=1 
  sh0_1 = ifelse(is.null(hyper$a0), 1, hyper$a0)
  sh0_0 = ifelse(is.null(hyper$b0), 1, hyper$b0)
  
  a = ifelse(is.null(hyper$c1), 1, hyper$c1)
  b = ifelse(is.null(hyper$c2), 1, hyper$c2)
  
  r = ifelse(is.null(hyper$d2), 1, hyper$d2)
  hatbeta=matrix(hat.r,ncol=p) + 10^-5
  progress = ifelse(debugging, 10^(floor(log10(iterations))-1), 0)
  if(robust){
  if(sparse){fit=BRGL_SS(xx1, y, CLC, p, d, iterations, hatAlpha, hatbeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq,
              xi1, xi2, r, a, b, sh0_1, sh0_0, progress)}
  else{fit=BRGL(xx1, y, CLC, p, d, iterations, hatAlpha, hatbeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq,
                xi1, xi2, r, a, b, progress)}
  }else{
    if(sparse){fit=BGLPointMass(xx1, y, CLC, p, d, iterations, hatAlpha, hat.r, invTAUsq.star, invSigAlpha0, hat.pi.s,
                                lambda.star, hat.sigma.sq, a.star, b.star, alpha, gamma, mu.star, nu.star, progress)}
    else{fit=BGL(xx1, y, CLC, p, d, iterations, hat.r, hatAlpha, invTAUsq.star, invSigAlpha0, lambda.star, 
                  hat.sigma.sq, a.star, b.star, alpha, gamma, progress)}
  }
  
}
