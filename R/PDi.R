
#' Phylogenetic Diversity indicator calculating function
#'
#' @inheritParams calculate_EDGE1
#' @inheritParams calculate_EDGE1_multiphylo
#' @inheritParams calculate_EDGE2
#' @inheritParams calculate_EDGE2_multiphylo
#' @inheritParams calculate_EDGE2_multiple
#'
#' @param RLcat.cols names of columns storing Red List assessments. to calculate Phylogenetic Diversity indicator.
#'  More than one column (generally used for different time points) are enabled, although PDi can be calculated for just one column.
#' @param ... More parameters to include in the EDGE2 calculating function (`ext.prob`, `summarise`, `seed`, `verbose`)
#'
#' @returns a data frame including where PD, ePDloss and PDi statistics.
#' If `summarise = FALSE`, a list with individual iteration and tree values will be returned.
#'
#'
#' @examplesIf interactive()
#'
#' PDi <-  calculate_PD_indicator(cycad.tree,
#'                                cycad.table,
#'                                RLcat.cols = c("RL_2003", "RL_2014"),
#'                                verbose = F,
#'                                summarise = T,
#'                                ext.prob ="IUCN50",
#'                                seed = 1234)
#'
#'
#' @author I. Ramos-Gutiérrez, R. Gumbs
#'
#' @export
#'
calculate_PD_indicator <- function(tree,
                                   table,
                                   species.col = "species",
                                   RLcat.cols,
                                   n.iter = NULL,
                                   summarise = TRUE,
                                   ...){

  if(!all(RLcat.cols %in% colnames(table))){stop("Please make sure all 'RLcat.cols' values are column names in 'table'")}

  match.arg(class(tree), c("phylo", "multiPhylo"))

  PDi_list <- vector(mode = "list", length = length(RLcat.cols))
  names(PDi_list) <- RLcat.cols

  for(time.col in RLcat.cols){
    table_tn <- table[,c(species.col, time.col)]
    colnames(table_tn) <- c("species", "RLcat")

    if(class(tree) == "phylo"){
      tree <- list(tree)
      class(tree) <- "multiPhylo"
    }

      edge_values_tn <- calculate_EDGE(tree,
                                       method = "EDGE2",
                                       table_tn,
                                       return.all = TRUE,
                                       summarise = FALSE,
                                       n.iter = n.iter,
                                       ...
      )




      if(length(tree) == 1 & is.null(n.iter)){
        edge_values_tn_comp <- edge_values_tn[[1]][["ePDloss"]]
      }else if(length(tree) == 1 & !is.null(n.iter)){
        edge_values_tn_comp <- lapply(1:n.iter, function(i){
          edge_values_tn[[i]][[1]][["ePDloss"]] |> dplyr::mutate(iter = i)
        })
      }else{
        edge_values_tn_comp <- lapply(1:n.iter, function(i){
          lapply(1:length(tree), function(t){
            edge_values_tn[[i]][[t]][["ePDloss"]] |> dplyr::mutate(iter = i, tree = t)
          })
        })
      }

      PDi_df <- edge_values_tn_comp |>
          dplyr::bind_rows() |>
          dplyr::mutate(dplyr::across(dplyr::everything(), unlist)) |>
          dplyr::mutate("PDi" = ePDloss / PD *100)


      if(isTRUE(summarise)){
        PDi_df <- PDi_df |>
          dplyr::summarise(PD_mn = mean(PD),
                           PD_sd = sd(PD),

                           PD_med = median(PD),
                           PD_iqr = IQR(PD),

                           ePDloss_mn = mean(ePDloss),
                           ePDloss_sd = sd(ePDloss),

                           ePDloss_med = median(ePDloss),
                           ePDloss_iqr = IQR(ePDloss),

                           PDi_mn = mean(PDi),
                           PDi_sd = sd(PDi),

                           PDi_med = median(PDi),
                           PDi_iqr = IQR(PDi),)
      }


    PDi_list[[time.col]] <- PDi_df
  }


  if(isTRUE(summarise)){
    PDi_list_sum <- lapply(RLcat.cols, function(RLcat){
    PDi_list[[RLcat]] |> dplyr::mutate(time.col = RLcat)
      }) |>
      bind_rows() |>
      relocate(time.col, .before = PD_mn)
    return(PDi_list_sum)
  }else{
    return(PDi_list)
  }
}

