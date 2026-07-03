#' Calculate EDGE1 index
#'
#' @param tree a phylo object.
#' @param table a tibble or data frame with two columns, including species names and Red List categories.
#' @param species.col column storing species names. If not specified, a column named "species" will be searched for.
#' @param RLcat.col column storing Red List assessments names. If not specified, a column named "RLcat" will be searched for.
#' @param sort.list Logical. If TRUE, the EDGE list will be sorted from higher to lower values.
#'
#' @return a data frame with the input table and the EDGE metric (sensu  Isaac et al., 2007).
#' @references Isaac, N.J., Turvey, S.T., Collen, B., Waterman, C. & Baillie, J.E. (2007)
#' Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2: e296.
#' \url{https://doi.org/10.1371/journal.pone.0000296}
#' @export
#'
calculate_EDGE1 <- function(tree,
                            table,
                            species.col = "species",
                            RLcat.col = "RLcat",
                            sort.list = TRUE,
                            ...){

  # check column names and rename
  if(!species.col %in% colnames(table)){
    stop("Column '", species.col, "' is not a column name in your table.\nPlease check or alternatively assign 'species' as column name in your table")
  }

  if(!RLcat.col %in% colnames(table)){
    stop("Column '", RLcat.col, "' is not a column name in your table.\nPlease check or alternatively assign 'RLcat' as column name in your table")
  }
  table <- table[,c(species.col, RLcat.col)]
  colnames(table) <-  c("species", "RLcat")

  if(!inherits(tree, "phylo")){stop("'tree' parameter must be an object of class 'phylo'")}
  if(!all(tree$tip.label %in% table$species)){
    warning("Some species in 'tree$tip.label' are not included in 'table$species'")
  }

  if(!all(table$species %in% tree$tip.label )){
    stop("Some species in 'table$species' are not included in 'tree$tip.label'")
  }


  if(!all(table$RLcat %in% c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW") )){
    stop("Categories should be: ", paste0(c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW"), collapse = " "))
  }

  table$pext <- NA_integer_
  table$pext[table$RLcat == "LC"] <- 0
  table$pext[table$RLcat == "NT"] <- 1
  table$pext[table$RLcat == "CD"] <- 1
  table$pext[table$RLcat == "VU"] <- 2
  table$pext[table$RLcat == "EN"] <- 3
  table$pext[table$RLcat == "CR"] <- 4
  table$pext[table$RLcat == "EW"] <- 4

  ED <- suppressWarnings(caper::ed.calc(tree))
  ED_res <- ED$spp

  table <- merge(table, ED_res)

  table$EDGE <- log(1+table$ED) + table$pext*(log(2))

  EDmed <- median(table$ED)
  table$isEDGEsp <- 0
  table$isEDGEsp[table$ED >= EDmed & table$RLcat %in% c("VU", "EN", "CR", "EW", "EX")] <- 1


  table <- table[,-which(colnames(table) == "pext")]

   if(isTRUE(sort.list)){
     table <- table[order(table$EDGE, decreasing = T),]
     rownames(table) <- 1:nrow(table)
   }


   return(table)
}


#' Calculate EDGE1 index for a multiphylo object
#'
#' @inheritParams calculate_EDGE1
#'
#' @param multitree a `multiPhylo` object storing multiple trees
#' @param summarise Logical. If FALSE, a list of dataframes outputted from `calculate_EDGE1` function will be returned,
#' one for each tree stored in the 'multitree' `multiPhylo` object. If FALSE, results will be summarised into a single table,
#' with mean, standard deviation, median and inter-quartile ranges for ED and EDGE metrics across trees.
#' @param parallelize Logical. If TRUE, several CPU cores (defined in `n.cores`) will be used.
#' @param n.cores Integer. Number of cores to use simultaneoulsly. Default is the number of available ones minus one.
#'
#' @return If `summarise` parameter is set to TRUE (default), a data frame with ED and EDGE averaged values. If not summarised,
#' a list storing results of `calculate_EDGE1` for each tree in the 'multitree' object.
#'
#' @references Isaac, N.J., Turvey, S.T., Collen, B., Waterman, C. & Baillie, J.E. (2007)
#' Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE 2: e296.
#' \url{https://doi.org/10.1371/journal.pone.0000296}
#' @export
#'
calculate_EDGE1_multiphylo <- function(multitree,
                                       table,
                                       species.col = "species",
                                       RLcat.col = "RLcat",
                                       sort.list = TRUE,
                                       #parameters for multi
                                       summarise = TRUE,
                                       parallelize = FALSE,
                                       n.cores = NULL,
                                       ...){

  # check column names and rename
  if(!species.col %in% colnames(table)){
    stop("Column '", species.col, "' is not a column name in your table.\nPlease check or alternatively assign 'species' as column name in your table")
  }

  if(!RLcat.col %in% colnames(table)){
    stop("Column '", RLcat.col, "' is not a column name in your table.\nPlease check or alternatively assign 'RLcat' as column name in your table")
  }
  table <- table[,c(species.col, RLcat.col)]
  colnames(table) <-  c("species", "RLcat")

  if(!inherits(multitree, "multiPhylo")){stop("'multitree' parameter must be an object of class 'multiPhylo'")}

  if(!all(unique(unlist(sapply(multitree, '[', "tip.label"))) %in% table$species)){
    warning("Some species in a 'tree$tip.label' are not included in 'table$species'")
  }

  if(!all(table$species %in% unique(unlist(sapply(multitree, '[', "tip.label"))))){
    stop("Some species in 'table$species' are not included in a 'tree$tip.label'")
  }


  if(!all(table$RLcat %in% c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW") )){
    stop("Categories should be: ", paste0(c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW"), collapse = " "))
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

    EDGElist <- future.apply::future_lapply(multitree,
                                            calculate_EDGE1,
                                            table = table,
                                            species.col = "species",
                                            RLcat.col = "RLcat",
                                            sort.list = sort.list,
                                            ...
    )
    future::plan(future::sequential)
  }else{
    EDGElist <- future.apply::future_lapply(multitree,
                                            calculate_EDGE1,
                                            table = table,
                                            species.col = "species",
                                            RLcat.col = "RLcat",
                                            sort.list = sort.list,
                                            ...
    )
  }

  EDGElist_tree <- lapply(1:length(multitree), function(i){EDGElist[[i]] |> dplyr::mutate(tree = i)})

  if (isTRUE(summarise)){

    EDGElist_ret <- EDGElist_tree |>
      dplyr::bind_rows() |>
      # dplyr::select(-`RLcat`) |>
      dplyr::group_by(species) |>
      dplyr::summarise(EDmed = median(ED),
                       EDiqr = IQR(ED),

                       EDGEmed = median(EDGE),
                       EDGEiqr = IQR(EDGE),

                       isEDGEsp = mean(isEDGEsp)) |>
      dplyr::mutate(isEDGEsp = ifelse(isEDGEsp >= 0.5, 1, 0)) |> # Over median in more than 50% of trees

      dplyr::left_join(table, by= "species") |>
      dplyr::relocate(RLcat , .after = species)

    if(isTRUE(sort.list)){EDGElist_ret <- dplyr::arrange(EDGElist_ret, dplyr::desc(EDGEmed))}

  }else {
    EDGElist_ret <- EDGElist_tree
  }


  return(EDGElist_ret)
}
