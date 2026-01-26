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



#' Multiple iteration EDGE calculating function.
#'
#' @param n.iter Integer. Number of times the function will be run.
#' @param parallelize Logical. If TRUE, several CPU cores will be used to compute EDGE scores.
#' @param n.cores Integer. Number of cores to use simultaneoulsly. Default is the number of available ones minus one.
#' @param seed Integer. seed number used to randomize iterations. if NULL, a random number will be enerated and printed
#' @param summarise Logical. IF TRUE (default), median and IQR values across iterations will be summarise. Only applicable when return.all is FALSE.
#'
#' @inheritParams calculate_EDGE2
#'
#' @returns A list of length = n.iter
#'
#' @author I. Ramos-Gutiérrez.
#'
#' @export
#'
calculate_EDGE_multiple <- function(tree,
                                     table,
                                     verbose = T,
                                     ext.prob = "Isaac",
                                     method ="EDGE2",
                                     sort.list = FALSE,
                                     return.all = FALSE,
                                     summarise = TRUE,
                                     n.iter = 10,
                                     parallelize = FALSE,
                                     n.cores = NULL,
                                     seed = NULL){


  if(is.null(seed)){
    seed <- round(runif(1, 1, 999999999))
    message(paste0("Seed has been set to: ", seed))
    set.seed(seed)

  }

  if(isTRUE(parallelize)){
    if(is.null(n.cores)){
      n.cores <- future::availableCores()-1
    }
    if(isTRUE(n.cores > future::availableCores())){
      message(paste0("n.cores value greater than available. Setting maximum-1 (", future::availableCores()-1, ")"))
      n.cores <- future::availableCores()-1
    }


    future::plan(future::multisession, workers = n.cores)

    EDGElist <- future.apply::future_lapply(1:n.iter,
                                            FUN = calculate_EDGE,
                                            method = method,
                                            tree = tree,
                                            table = table,
                                            ext.prob = ext.prob,
                                            sort.list = sort.list,
                                            return.all = return.all,

                                            future.seed = seed)
    future::plan(future::sequential)
  }else{
    EDGElist <- future.apply::future_lapply(1:n.iter,
                                            FUN = calculate_EDGE,
                                            method = method,
                                            tree = tree,
                                            table = table,
                                            ext.prob = ext.prob,
                                            sort.list = sort.list,
                                            return.all = return.all,

                                            future.seed = seed)
  }

  if (isFALSE(return.all) & isTRUE(summarise)){

    EDGElist_compl <- dplyr::bind_rows(EDGElist) |>
      dplyr::group_by(species) |>
      dplyr::summarise(RL.cat = unique(RL.cat),

                       TBLmed = unique(TBL),

                       pextmed = median(pext),
                       pextiqr = IQR(pext),

                       EDmed = median(ED),
                       EDiqr = IQR(ED),

                       EDGEmed = median(EDGE),
                       EDGEiqr = IQR(EDGE)
      )
  }else{
    EDGElist_compl <- EDGElist
  }
  return(EDGElist_compl)

}






#' EDGE2 calculating function for a multiPhylo tree object.
#'
#'
#' @inheritParams calculate_EDGE_multiple
#'
#' @returns A list of length equal to the number of trees included in the multiphylo object, storing
#' information resulting from the calculate_EDGE function used based on the specified method.
#'
#' @author I. Ramos-Gutiérrez.
#'
#' @export
#'
calculate_EDGE_multiphylo <- function(multiphylo,
                                      table,
                                      method = "EDGE2",
                                      ext.prob = "Isaac",
                                      sort.list = FALSE,
                                      return.all = FALSE,
                                      summarise = TRUE,
                                      parallelize = FALSE,
                                      n.cores = NULL,
                                      seed = NULL,
                                      verbose = TRUE){

  if(!inherits(multiphylo, "multiPhylo")){stop("multiphylo should be an obhecto of class 'multiPhylo'")}
  if(is.null(seed)){
    seed <- round(runif(1, 1, 999999999))
    message(paste0("Seed has been set to: ", seed))
    set.seed(seed)

  }
  # parallelize <- FALSE # STILL UNDER CONSTRUCTION; the ... are being troubly
  if(isTRUE(parallelize)){
    if(is.null(n.cores)){
      n.cores <- future::availableCores()-1
    }
    if(isTRUE(n.cores > future::availableCores())){
      message(paste0("n.cores value greater than available. Setting maximum-1 (", future::availableCores()-1, ")"))
      n.cores <- future::availableCores()-1
    }


    future::plan(future::multisession, workers = n.cores)

    EDGElist <- future.apply::future_lapply(multiphylo,
                                            calculate_EDGE,
                                            table = table,
                                            method = method,
                                            ext.prob = ext.prob,
                                            sort.list = sort.list,
                                            return.all = return.all,
                                            verbose = FALSE,

                                            future.seed = seed
    )
    future::plan(future::sequential)
    }else{
    EDGElist <- future.apply::future_lapply(multiphylo,
                                            calculate_EDGE,
                                            table = table,
                                            method = method,
                                            ext.prob = ext.prob,
                                            sort.list = sort.list,
                                            return.all = return.all,
                                            verbose = verbose,

                                                 future.seed = seed
    )
  }

  if (isFALSE(return.all) & isTRUE(summarise)){

    EDGElist_compl <- dplyr::bind_rows(EDGElist) |>
      dplyr::group_by(species) |>
      dplyr::summarise(RL.cat = unique(RL.cat),

                       TBLmed = unique(TBL),

                       pextmed = median(pext),
                       pextiqr = IQR(pext),

                       EDmed = median(ED),
                       EDiqr = IQR(ED),

                       EDGEmed = median(EDGE),
                       EDGEiqr = IQR(EDGE)
      )
  }else{
    EDGElist_compl <- EDGElist
  }
  return(EDGElist_compl)

}
