#' Variable selection for a pqrBayes object
#'
#' Variable selection for a pqrBayes object
#'
#' @param object a pqrBayes object.
#' @param sparse logical flag. If TRUE, the sparse model is used for variable selection. The default value is TRUE.
#' @param model the model to be fitted. Users can also choose "linear" for a linear model or "VC" for a varying coefficient model.
#'
#' @details For class `Sparse', the median probability model (MPM) (Barbieri and Berger, 2004) is used to identify predictors that are significantly associated
#' with the response variable. For class `NonSparse', variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2023). Robust Bayesian variable selection for gene-environment interactions. 
#' {\emph{Biometrics}, 79(2), 684-694} \doi{10.1111/biom.13670}
#'
#' Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model selection. {\emph{Ann. Statist}, 32(3):870–897}
#'
#' @rdname select.pqrBayes
#' @usage pqrBayes.select(object,sparse=T,model="linear")
#' @return an object of class `select' is returned, which includes the indices of the selected predictors (e.g. genetic factors).
#'
#' @seealso \code{\link{pqrBayes}}
#'
#' @examples
#' ## The quantile regression model
#' data(data)
#' data = data$data_linear
#' g=data$g
#' y=data$y
#' e=data$e
#' 
#' fit1=pqrBayes(g,y,u=NULL,e,quant=0.5,kn=NULL,degree=NULL,model="linear")
#' sparse=TRUE
#' select=pqrBayes.select(obj = fit1,sparse = sparse,model="linear")
#' 
#' ## The quantile varying coefficient model
#' data(data)
#' data = data$data_varying
#' g=data$g
#' y=data$y
#' u=data$u
#' e=data$e
#' fit1=pqrBayes(g,y,u,e,kn=2,degree=2,quant=0.5,model="VC")
#' sparse=TRUE
#' select=pqrBayes.select(obj = fit1,sparse = sparse,model="VC")
#' select
#'
#' \donttest{
#' ## non-sparse
#' sparse=FALSE
#' fit2=pqrBayes(g,y,u,e,quant=0.5,kn=2,degree=2,sparse = sparse,model="VC")
#' select=pqrBayes.select(obj=fit2,sparse=FALSE,model="VC")
#' select
#' }
#'
#' @export
pqrBayes.select <- function(object,sparse=T,model="linear"){
  if(model=="VC"){
    select = VCselect(obj = object,sparse = sparse)
  }
  else if(model=="linear"){
    select = linselect(obj = object,sparse = sparse)
  }
  else{
    stop("model should be either VC or linear")
  }
  class(select)="select"
  return(select)
}
  
  
 