#' Calculate EDGE wrapping function

#' @param method Method to be used. It has to be either "EDGE1" or "EDGE2" (default).
#' @param tree Phylogenetic information as a 'phylo' or 'multiPhylo' object
#' @inheritParams calculate_EDGE1
#' @inheritParams calculate_EDGE2
#' @inheritParams calculate_EDGE1_multiphylo
#' @inheritParams calculate_EDGE2_multiphylo
#'
#'
#' @examplesIf interactive()
#'
#' # Calculate EDGE-1 for a single tree
#' calculate_EDGE(monotreme.tree,
#' monotreme.table,
#' method = "EDGE1")
#'
#' # Calculate EDGE-1 for a multiPhylo object (multiple phylogenies)
#' calculate_EDGE(cycad.multitree, cycad.table,
#'                method = "EDGE1",
#'                RLcat.col = "RL_2003", # in this example it is not the default name (RLcat)!
#'                sort.list = T)
#'
#' # Other variations
#' calculate_EDGE(cycad.multitree,
#'                cycad.table,
#'                method = "EDGE1",
#'                RLcat.col = "RL_2014",
#'                sort.list = FALSE,
#'                summarise = FALSE,
#'                parallelize = TRUE,
#'                n.cores = 4)
#'
#'
#' # Calculate EDGE-2 for a single tree
#' calculate_EDGE(monotreme.tree,
#'                monotreme.table,
#'                method = "EDGE2",
#'                seed=123)
#'
#' calculate_EDGE(monotreme.tree,
#'                monotreme.table,
#'                method = "EDGE2",
#'                sort.list = FALSE,
#'                ext.prob = "IUCN100",
#'                seed=123)
#'
#' # Calculate EDGE-2 for a multiPhylo object (multiple phylogenies)
#' calculate_EDGE(cycad.multitree,
#'                cycad.table,
#'                method = "EDGE2",
#'                RLcat.col = "RL_2014",
#'                sort.list = FALSE,
#'                summarise = FALSE,
#'                parallelize = TRUE,
#'                n.cores = 4)
#'
#' # Calculate EDGE-2 for a single tree, but 100 times
#' calculate_EDGE(monotreme.tree,
#'                monotreme.table,
#'                n.iter = 100)
#'
#' # Calculate EDGE-2 for a multiple phylogenies, and 10 times each, combining other parameters
#' calculate_EDGE(tree = cycad.multitree,
#'                table = cycad.table,
#'                RLcat.col = "RL_2003",
#'                ext.prob ="IUCN50",
#'                n.iter = 10,
#'                parallelize = T,
#'                sort.list = T,
#'                return.all = F,
#'                summarise = F)
#'
#' @export
#'
calculate_EDGE <- function(tree, table, method = "EDGE2", n.iter = NULL, ...){

  if(!method %in% c("EDGE1", "EDGE2")){
    stop("'method' argument has to be \"EDGE1\" or \"EDGE2\"")
  }
  if(!(class(tree) %in% c("phylo", "multiPhylo"))){
    stop("'tree' should be an object of class 'phylo' or 'multiPhylo'")
  }


  if(method == "EDGE1" & class(tree) == "phylo"){
    result <- calculate_EDGE1(tree, table, ...)
  }

  if(method == "EDGE1" & class(tree) == "multiPhylo"){
    result <- calculate_EDGE1_multiphylo(tree, table, ...)
  }

  if(method == "EDGE2" & class(tree) == "phylo" & is.null(n.iter)){
    result <- calculate_EDGE2(tree, table, ...)
  }

  if(method == "EDGE2" & class(tree) == "multiPhylo" & is.null(n.iter)){
    result <- calculate_EDGE2_multiphylo(tree, table, ...)
  }

  if(method == "EDGE2" & !is.null(n.iter)){
    result <- calculate_EDGE2_multiple(tree, table, n.iter = n.iter, ...)
  }

  return(result)
}


