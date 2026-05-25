
#' Phylogenetic Diversity indicator calculating function
#'
#' @param tree An object of class 'phylo' or 'mutliPhylo' where to calculate the PDi
#' @param table a data.frame containing a column names 'species' with all taxa in the tree,
#' and columns storing Red List categories.
#' @param time.cols column names of all times whose PDi should be calculated.
#'
#'
#' @returns a list including a dataframe for each time column, where PD, ePDloss and PDi are shown
#'
#'
#' @examples
#' calculate_PD_indicator(trees = cycad.tree,
#'                        table = cycad.table,
#'                        time.cols = c("RL_2003", "RL_2014"))
#'
#'
#' @author I. Ramos-Gutiérrez, R. Gumbs
#'
#' @export
#'
calculate_PD_indicator <- function(tree, table, time.cols, seed = NULL, ...){

  if(!all(time.cols %in% colnames(table))){stop("Please make sure all 'time.cols' values are column names in 'table'")}
  match.arg(class(tree), c("phylo", "multiPhylo"))

  PDi_list <- vector(mode = "list", length = length(time.cols))
  names(PDi_list) <- time.cols

  for(time.col in time.cols){
    table_tn <- table[,c("species", time.col)]
    colnames(table_tn) <- c("species", "RLcat")

    if(class(tree) == "phylo"){
      tree <- list(tree)
      class(tree) <- "multiPhylo"
    }

      edge_values_tn <- calculate_EDGE2_multiphylo(tree,
                                                  table_tn,
                                                  return.all = TRUE,
                                                  summarise = FALSE,
                                                  seed = seed,
                                                  ...
      )




    epdl_vals_tn <- sapply(edge_values_tn, "[[", "ePDloss") |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(dplyr::across(dplyr::everything(), unlist)) |>
      dplyr::mutate("PDi" = ePDloss / PD *100)


    PDi_list[[time.col]] <- epdl_vals_tn
  }

  return(PDi_list)
}

