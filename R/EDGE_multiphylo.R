
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
                                      sort.list = FALSE,
                                      # ext.prob = "Isaac",
                                      # return.all = FALSE,
                                      summarise = TRUE,
                                      parallelize = FALSE,
                                      n.cores = NULL,
                                      seed = NULL,
                                      # verbose = TRUE,
                                      ...){

  if(!inherits(multiphylo, "multiPhylo")){stop("multiphylo should be an object of class 'multiPhylo'")}
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

    EDGElist <- future.apply::future_lapply(multiphylo,
                                            calculate_EDGE,
                                            table = table,
                                            method = method,
                                            sort.list = sort.list,
                                            # ext.prob = ext.prob,      #just for EDGE2, goes to dots
                                            # return.all = return.all,  #just for EDGE2, goes to dots
                                            # verbose = FALSE,
                                            ... ,

                                            future.seed = seed
    )
    future::plan(future::sequential)
  }else{
    EDGElist <- future.apply::future_lapply(multiphylo,
                                            calculate_EDGE,
                                            table = table,
                                            method = method,
                                            sort.list = sort.list,
                                            # ext.prob = ext.prob,      #just for EDGE2, goes to dots
                                            # return.all = return.all,  #just for EDGE2, goes to dots
                                            # verbose = FALSE,
                                            ... ,

                                            future.seed = seed
    )
  }

  if(!exists("return.all")){return.all <- FALSE}

  if (!isTRUE(return.all) & isTRUE(summarise)){

    EDGElist_compl <- dplyr::bind_rows(EDGElist) |>
      dplyr::select(-`RL.cat`) |>
      dplyr::group_by(species) |>
      dplyr::summarise(TBLmed = mean(TBL),

                       pextmed = median(pext),
                       pextiqr = IQR(pext),

                       EDmed = median(ED),
                       EDiqr = IQR(ED),

                       EDGEmed = median(EDGE),
                       EDGEiqr = IQR(EDGE)
      ) |>
    dplyr::left_join(table, by= "species") |>
    dplyr::relocate(RL.cat , .after = species)

    if(isTRUE(sort.list)){EDGElist_compl <- dplyr::arrange(EDGElist_compl, dplyr::desc(EDGEmed))}
  }else{
    EDGElist_compl <- EDGElist
  }
  return(EDGElist_compl)

}
