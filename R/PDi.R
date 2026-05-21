
#' Phylogenetic Diversity indicator calculating function
#'
#' @param trees An object of class 'phylo' or 'mutliPhylo' where to calculate the PDi
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
calculate_PD_indicator <- function(trees, table, time.cols, ...){

  if(!all(time.cols %in% colnames(table))){stop("Please make sure all 'time.cols' values are column names in 'table'")}
  match.arg(class(trees), c("phylo", "multiPhylo"))

  PDi_list <- vector(mode = "list", length = length(time.cols))
  names(PDi_list) <- time.cols

  for(time.col in time.cols){
    table_tn <- table[,c("species", time.col)]
    colnames(table_tn) <- c("species", "RL.cat")

    if(class(trees) == "phylo"){
      trees <- list(trees)
      class(trees) <- "multiPhylo"
    }

      edge_values_tn <- calculate_EDGE_multiphylo(trees,
                                                  table_tn,
                                                  return.all = TRUE,
                                                  summarise = FALSE,
                                                  ...
      )




    epdl_vals_tn <- sapply(edge_values_tn, "[[", 3) |>
      t() |>
      as.data.frame() |>
      dplyr::mutate(dplyr::across(dplyr::everything(), unlist)) |>
      dplyr::mutate("PDi" = ePDloss / PD *100)


    PDi_list[[time.col]] <- epdl_vals_tn
  }

  return(PDi_list)
}

