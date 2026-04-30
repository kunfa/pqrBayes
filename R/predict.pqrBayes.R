#' Make predictions from a pqrBayes object
#'
#' Make predictions from a pqrBayes object
#'
#' @param object a pqrBayes object.
#' @param g.new a matrix of new predictors (e.g. genetic factors) at which predictions are to be made. When being applied to the sparse linear, binary LASSO or group LASSO, g.new = g.
#' @param u.new a vector of new environmental factor at which predictions are to be made. When being applied to the sparse linear model, binary LASSO or group LASSO, u.new = NULL. The default value is NULL.
#' @param e.new a vector or matrix of new clinical covariates at which predictions are to be made. The default value is NULL.
#' @param y.new a vector of the response of new observations. When being applied to the sparse linear model, binary LASSO or group LASSO, y.new = y.
#' @param robust logical flag. If TRUE, robust methods are used. Otherwise, non-robust methods are used. The default value is TRUE.
#' @param quant the quantile level specified by users. Required when robust = TRUE. Ignored (set to NULL) when robust = FALSE.The default value is 0.5.
#' @param model the model to be fitted. The default is "VC" for a quantile varying coefficient model. Users can also specify "linear" for a sparse linear model, "binary" for binary LASSO and "group" for a group LASSO.
#' @param ... other predict arguments
#' 
#' @details g.new (u.new) must have the same number of columns as g (u) used for fitting the model. By default, the clinical covariates are NULL unless 
#' provided. The predictions are made based on the posterior estimates of coefficients in the pqrBayes object.
#'
#'
#' @usage predict_pqrBayes(object, g.new, u.new, e.new, y.new, robust, quant, model, ...)
#' @return  an object of class `pqrBayes.pred' is returned, which is a list with components:
#' \item{error}{prediction error.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict_pqrBayes
#' @seealso \code{\link{pqrBayes}}
#' @examples
#' ## The quantile regression model
#' data(data)
#' data = data$data_linear
#' g=data$g
#' y=data$y
#' e=data$e
#' fit1=pqrBayes(g,y,e,model="linear")
#' prediction=predict_pqrBayes(fit1,g, y.new = y, model="linear")
#' @export
predict_pqrBayes=function(object, g.new, u.new=NULL, e.new=NULL, y.new, robust=T, quant=0.5,model,...){
  if (model != "VC" && !is.null(u.new)) {
    
    stop("Argument 'u.new' should be NULL unless model = 'VC'.")
    
  }
  
  if (model == "VC" && is.null(u.new)) {
    
    stop("Argument 'u.new' must be provided when model = 'VC'.")
    
  }
  # quant check
  if (robust) {
    if (is.null(quant)) {
      stop("quant must be specified when robust = TRUE.")
    }
  } else {
    if (!is.null(quant)) {
      stop("quant must be NULL when robust = FALSE.")
    }
  }
  if(model=="VC"){
    pqrBayes.pred = predict_vc(object, g.new, u.new, e.new, y.new, robust, quant,...)
  }else if(model=="linear"){
    pqrBayes.pred = predict_lin(object, g.new, e.new, y.new, robust, quant,...)
  }else if(model=="binary"){
    pqrBayes.pred = predict_bin(object, g.new, e.new, y.new,...)
  }else if(model=="group"){
    pqrBayes.pred = predict_lin(object, g.new, e.new, y.new, robust, quant,...)
  }
  
  else{
    stop("model should be either VC, linear, binary or group")
  }
  class(pqrBayes.pred) = "pqrBayes.pred"
  return(pqrBayes.pred)
  #pred
}

