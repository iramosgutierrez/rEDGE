#' Calculate EDGE wrapping function

#' @param method Method to be used. It has to be either "EDGE1" or "EDGE2".
#' The output will be the result of running calculate_EDGE1, or calculate_EDGE2 depending on the method specified.
#' @inheritParams calculate_EDGE2
#'
#' @export
#'
calculate_EDGE <- function(tree, table, method = "EDGE2", ...){

  if(!method %in% c("EDGE1", "EDGE2")){
    stop("'method' argument has to be \"EDGE1\" or \"EDGE2\"")
  }

  if(method == "EDGE1"){
    result <- calculate_EDGE1(tree, table, ...)
  }else if(method == "EDGE2"){
    result <- calculate_EDGE2(tree, table, ...)
  }
  return(result)
}

