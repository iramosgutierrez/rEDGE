

#' Multiple iteration EDGE2 calculating function.
#'
#' @inheritParams calculate_EDGE2
#' @inheritParams calculate_EDGE2_multiphylo
#'
#' @param n.iter Integer. Number of times the function will be run.
#' @param tree A 'phylo' or 'multiPhylo' object.
#'
#'
#' @returns A list of length = n.iter
#'
#' @examplesIf interactive()
#' calculate_EDGE2_multiple(monotreme.tree,
#'                         monotreme.table,
#'                         verbose = T,
#'                         sort.list = T,
#'                         n.iter = 4,
#'                         parallelize = T,
#'                         n.cores = 4,
#'                         seed = 123456)
#' # example code
#'
#' @author I. Ramos-Gutiérrez.
#'
#' @export
#'
calculate_EDGE2_multiple <- function(tree,
                                     table,
                                     species.col = "species",
                                     RLcat.col = "RLcat",
                                     sort.list = FALSE,
                                     ext.prob = "Isaac",
                                     return.all = FALSE,
                                     summarise = TRUE,
                                     n.iter = 10,
                                     parallelize = FALSE,
                                     n.cores = NULL,
                                     seed = NULL,
                                     verbose = TRUE){


  if(!(class(tree) %in% c("phylo", "multiPhylo"))){stop("'tree' should be an object of class 'phylo' or 'multiPhylo'")}

  if(is.null(seed)){
    seed <- round(runif(1, 1, 999999999))
    if(isTRUE(verbose)){message(paste0("Seed for `future` package has been set to: ", seed, "\n"))}
    set.seed(seed)
  }

  table <- table[,c(species.col, RLcat.col)]
  colnames(table) <-  c("species", "RLcat")


  if(isTRUE(parallelize)){
    if(is.null(n.cores)){
      n.cores <- future::availableCores()-1
    }
    if(isTRUE(n.cores > future::availableCores())){
      message(paste0("n.cores value greater than available. Setting maximum-1 (", future::availableCores()-1, ")"))
      n.cores <- future::availableCores()-1
    }
  }


    if(isTRUE(parallelize)){
    future::plan(future::multisession, workers = n.cores)
    }

    if(inherits(tree, "phylo")){
    EDGElist <- future.apply::future_lapply(1:n.iter,
                                            FUN = calculate_EDGE2,
                                            tree = tree,
                                            table = table,
                                            species.col = "species",
                                            RLcat.col = "RLcat",
                                            ext.prob = ext.prob,
                                            sort.list = sort.list,
                                            return.all = return.all,
                                            verbose = verbose,
                                            # ...,

                                            future.seed = seed)
    }else{
      EDGElist <- future.apply::future_lapply(1:n.iter,
                                              FUN = calculate_EDGE2_multiphylo,
                                              multitree = tree,
                                              table = table,
                                              species.col = "species",
                                              RLcat.col = "RLcat",
                                              ext.prob = ext.prob,
                                              sort.list = sort.list,
                                              return.all = return.all,
                                              verbose = verbose,
                                              summarise = FALSE, # In any case, summarise after!
                                              parallelize = FALSE, # To avoid double parallellization!
                                              # ...,

                                              future.seed = seed)
    }

    if(isTRUE(parallelize)){
    future::plan(future::sequential)
    }

if(!isTRUE(return.all)){
    if(inherits(tree, "phylo")){
      EDGElist_tree <- lapply(1:n.iter, function(i){
          EDGElist[[i]] |> dplyr::mutate(iter = i)
        })
    }else{
    EDGElist_tree <- lapply(1:n.iter, function(i){
                                  lapply(1:length(tree), function(t){
                                    EDGElist[[i]][[t]] |> dplyr::mutate(iter = i)
                                  })
                 })
    }
}

  if(!isTRUE(return.all) & isTRUE(summarise)){

    EDGElist_ret <- EDGElist_tree |>
      dplyr::bind_rows() |>
      dplyr::group_by(species) |>
      dplyr::summarise(TBLmn = mean(TBL),

                       pextmed = median(pext),
                       pextiqr = IQR(pext),
                       pextmed = median(pext),
                       pextiqr = IQR(pext),

                       EDmn = mean(ED),
                       EDsd = sd(ED),
                       EDmed = median(ED),
                       EDiqr = IQR(ED),

                       EDGEmn = mean(EDGE),
                       EDGEsd = sd(EDGE),
                       EDGEmed = median(EDGE),
                       EDGEiqr = IQR(EDGE)     ) |>
      dplyr::left_join(table, by= "species") |>
      dplyr::relocate(RLcat , .after = species)

    if(isTRUE(sort.list)){EDGElist_ret <- dplyr::arrange(EDGElist_ret, dplyr::desc(EDGEmed))}
}
    if( isTRUE(return.all)                     ){return(EDGElist)}
    if(!isTRUE(return.all) & !isTRUE(summarise)){return(EDGElist_tree)}
    if(!isTRUE(return.all) &  isTRUE(summarise)){return(EDGElist_ret)}

}



