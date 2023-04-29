#' #' @rdname refinedmod
#' #' @export
#' #' @param object an object class of \code{refinedmod}
#' #' @param model index indicating the summary call to print
#' #' @param ... other argument not in use at the moment 
#' summary.refinedmod <- function(object, model = NULL, ...){
#'   if(!is.null(model)){
#'     if(model==1){
#'       out <- capture.output(object[["allsummary.coef"]][["out.direct"]])
#'       v <- c("Model summary of the direct estimation:", out, "\n") # char vec
#'       s <- paste(v, collapse = "\n") # single string
#'     }else if(model == 2){
#'       out <- capture.output(object[["allsummary.coef"]][["out.inverse"]])
#'       v <- c("Model summary of the inverse estimation:", out, "\n") # char vec
#'       s <- paste(v, collapse = "\n") # single string
#'     }else{
#'       out <- capture.output(object[["allsummary.coef"]])
#'       v <- c("Model summary of the direct, inverse and inverse estimation:", out, "\n") # char vec
#'       s <- paste(v, collapse = "\n")}
#'   }else{
#'     out <- capture.output(object[["allsummary.coef"]][["out.gee"]])
#'     v <- c("Model summary of the refined GEE estimation:", out, "\n") # char vec
#'     s <- paste(v, collapse = "\n") # single string
#'   }
#'   # cat(s)
#'   structure(s, class = "summary.refinedmod")
#' }
#' 
#' 
#' #' @export 
#' print.summary.refinedmod <- function(x,...){
#'   cat(x)
#'   invisible(x)
#' }
#' 
#' #' #' @export 
#' #'  print.summary2.refinedmod <- function(x,...){
#' #'    cat(attr(x, "message"))
#' #'    print.data.frame(x)
#' #'    invisible(x)
#' #'  }