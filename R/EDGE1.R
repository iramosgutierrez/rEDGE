#' Calculate EDGE index
#'
#' @param tree a phylo object.
#' @param table a tibble or data frame with two columns, named "species" and "RL.cat"
#'
#' @return a data frame with the input table and the EDGE metric (sensu  Isaac et al., 2007).
#' @export
#'
calculate_EDGE1 <- function(tree, table){

  if(!all(tree$tip.label %in% table$species )){
    stop("Some species in 'tree$tip.label' are not included in 'table$species'")
  }

  if(!all(table$species %in% tree$tip.label )){
    stop("Some species in 'table$species' are not included in 'tree$tip.label'")
  }


  if(!all(table$RL.cat %in% cat_pext()$rl.cat )){
    stop("Categories should be: ", paste0(cat_pext()$rl.cat, collapse = " "))
  }

  table$pext[table$RL.cat == "LC"] <- 0
  table$pext[table$RL.cat == "NT"] <- 1
  table$pext[table$RL.cat == "VU"] <- 2
  table$pext[table$RL.cat == "EN"] <- 3
  table$pext[table$RL.cat == "CR"] <- 4


  ED <- caper::ed.calc(tree)
  ED_res <- ED$spp

  table <- merge(table, ED_res)

  table$EDGE <- log(1+table$ED) + table$pext*(log(2))

   table <- table[,-which(colnames(table) == "pext")]

   return(table)
}
