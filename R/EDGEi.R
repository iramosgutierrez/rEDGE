
#' EDGE index calculating function
#'
#' @inheritParams calculate_EDGE1
#' @inheritParams calculate_EDGE1_multiphylo
#' @inheritParams calculate_EDGE2
#' @inheritParams calculate_EDGE2_multiphylo
#' @inheritParams calculate_EDGE2_multiple
#' @inheritParams calculate_PD_indicator
#'
#' @returns a data frame
#' If `summarise = FALSE`, a list with individual iteration and tree values will be returned.
#'
#'
#' @examplesIf interactive()
#'
#' PDi <-  calculate_EDGE_index(tree = cycad.tree,
#'                              table = cycad.table,
#'                              RLcat.cols = c("RL_2003", "RL_2014"),
#'                              verbose = F,
#'                              summarise = T,
#'                              ext.prob ="IUCN50",
#'                              seed = 1234)
#'
#' @author I. Ramos-Gutiérrez
#'
#' @export
calculate_EDGE_index <- function(tree,
                                 table,
                                 species.col = "species",
                                 RLcat.cols,
                                 n.iter = 10,
                                 summarise = TRUE,
                                 ...){

  if(!all(RLcat.cols %in% colnames(table))){stop("Please make sure all 'RLcat.cols' values are column names in 'table'")}

  match.arg(class(tree), c("phylo", "multiPhylo"))

  EDGEi_list <- vector(mode = "list", length = length(RLcat.cols))
  names(EDGEi_list) <- RLcat.cols

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
                                     return.all = FALSE,
                                     summarise = FALSE,
                                     n.iter = n.iter,
                                     ...
    )


    if(length(tree) == 1 & is.null(n.iter)){

      edge_spp_tn_comp <- edge_values_tn[[1]]|>
        dplyr::summarise(nspp = dplyr::n(),
                         nEDGEspp = sum(isEDGEsp == 1, na.rm = TRUE),
                         EDGEi = nEDGEspp / nspp)

    }else if(length(tree) == 1 & !is.null(n.iter)){

      edge_spp_tn_comp <- lapply(1:n.iter, function(i){
        edge_spp_tn_comp <- edge_values_tn[[i]][[1]]|>
          dplyr::summarise(nspp = dplyr::n(),
                           nEDGEspp = sum(isEDGEsp == 1, na.rm = TRUE),
                           EDGEi = nEDGEspp / nspp) |>
          dplyr::mutate(iter = i)
      })

    }else{

      edge_spp_tn_comp <- lapply(1:n.iter, function(i){
        lapply(1:length(tree), function(t){
          edge_values_tn[[i]][[t]] |>
            dplyr::summarise(nspp = dplyr::n(),
                             nEDGEspp = sum(isEDGEsp == 1, na.rm = TRUE),
                             EDGEi = nEDGEspp / nspp) |>
            dplyr::mutate(iter = i, tree = t)
        })
      })

    }


    EDGEi_df <- edge_spp_tn_comp |>
          dplyr::bind_rows() |>
          dplyr::mutate(dplyr::across(dplyr::everything(), unlist))


      if(isTRUE(summarise)){
        EDGEi_df <- EDGEi_df |>
          dplyr::summarise(nspp  = mean(nspp ),

                           nEDGEspp_mn = mean(nEDGEspp),
                           nEDGEspp_sd = sd(nEDGEspp),


                           EDGEi_mn = mean(EDGEi),
                           EDGEi_sd = sd(EDGEi)             )
      }


    EDGEi_list[[time.col]] <- EDGEi_df

  }



  if(isTRUE(summarise)){
    EDGEi_list_sum <- lapply(RLcat.cols, function(RLcat.col){
      EDGEi_list[[RLcat.col]] |> dplyr::mutate(time.col = RLcat.col)
    }) |>
      dplyr::bind_rows() |>
      dplyr::relocate(time.col, .before = nspp)
    return(EDGEi_list_sum)
  }else{
    return(EDGEi_list)
  }
}
