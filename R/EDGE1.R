#' Calculate EDGE index
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
                            sort.list = FALSE,
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

  table <- table[,-which(colnames(table) == "pext")]

   if(isTRUE(sort.list)){
     table <- table[order(table$EDGE, decreasing = T),]
     rownames(table) <- 1:nrow(table)
   }

   # EDmedian <- median(table$ED)
   # table$EDGEspp <- "N"
   # table$EDGEspp[table$ED > EDmedian & table$RL.cat %in% c("VU", "EN", "CR", "EW")]

   return(table)
}


calculate_EDGE1_multiphylo <- function(multitree,
                                       table,
                                       species.col = "species",
                                       RLcat.col = "RLcat",
                                       sort.list = FALSE,
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

  if(!all(unique(unlist(sapply(multitree, '[', 4))) %in% table$species)){
    warning("Some species in a 'tree$tip.label' are not included in 'table$species'")
  }

  if(!all(table$species %in% unique(unlist(sapply(multitree, '[', 4))))){
    stop("Some species in 'table$species' are not included in a 'tree$tip.label'")
  }


  if(!all(table$RLcat %in% c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW") )){
    stop("Categories should be: ", paste0(c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW"), collapse = " "))
  }


  lapply(multitree, calculate_EDGE1, table = table,
         species.col = species.col,
         RLcat.col = RLcat.col,
         sort.list = sort.list)
  # EDmedian <- median(table$ED)
  # table$EDGEspp <- "N"
  # table$EDGEspp[table$ED > EDmedian & table$RL.cat %in% c("VU", "EN", "CR", "EW")]

  return(table)
}
