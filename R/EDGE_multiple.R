

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



