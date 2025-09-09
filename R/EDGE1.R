#' Calculate EDGE index
#'
#' @param tree a phylo object.
#' @param table a tibble or data frame with two columns, named "species" and "RL.cat"
#' @param sort.list Logical. If TRUE, the EDGE list will be sorted from higher to lower values.
#'
#' @return a data frame with the input table and the EDGE metric (sensu  Isaac et al., 2007).
#' @export
#'
calculate_EDGE1 <- function(tree,
                            table,
                            sort.list = FALSE){

  if(!all(tree$tip.label %in% table$species )){
    stop("Some species in 'tree$tip.label' are not included in 'table$species'")
  }

  if(!all(table$species %in% tree$tip.label )){
    stop("Some species in 'table$species' are not included in 'tree$tip.label'")
  }


  if(!all(table$RL.cat %in% c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW") )){
    stop("Categories should be: ", paste0(c(cat_pext()$rl.cat, "CD", "NE", "DD", "EW"), collapse = " "))
  }

  table$pext <- NA_integer_
  table$pext[table$RL.cat == "LC"] <- 0
  table$pext[table$RL.cat == "NT"] <- 1
  table$pext[table$RL.cat == "CD"] <- 1
  table$pext[table$RL.cat == "VU"] <- 2
  table$pext[table$RL.cat == "EN"] <- 3
  table$pext[table$RL.cat == "CR"] <- 4
  table$pext[table$RL.cat == "EW"] <- 4

  ED <- suppressWarnings(caper::ed.calc(tree))
  ED_res <- ED$spp

  table <- merge(table, ED_res)

  table$EDGE <- log(1+table$ED) + table$pext*(log(2))

   table <- table[,-which(colnames(table) == "pext")]

   if(isTRUE(sort.list)){
     table <- table[order(table$EDGE, decreasing = T),]
   }

   return(table)
}
