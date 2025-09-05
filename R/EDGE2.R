# remove excess pext, and reorders to same order as tree$tip.label
into_order <- function(tree, pext){
  new_pext <- pext[match(tree$tip.label, pext$species),]
  return (new_pext)
}

# order tree components
reorder_tree <- function(tree, ordering){
  tree@edge.length <- tree@edge.length[ordering]
  tree@edge <- tree@edge[ordering,]
  return(tree)
}

#' EDGE2 calculating function.
#'
#' @param verbose Logical. Should progress be printed or not.
#'
#' @inheritParams calculate_EDGE1
#'
#' @returns A list containing 3 slots.
#' \itemize{
#'  \item{list: }{will output specific TBL (terminal branch length), pext (sampled probability of extinction), ED (evolutionary distinctiveness), and EDGE (EDGE index).}
#'  \item{tree: }{the input tree.}
#'  \item{ePDloss: }{complete PD (phylogenetic diversity) and ePDloss (expected PD loss). Units are as in the input tree  (generally Million years).}
#' }
#'
#' @author I. Ramos-Gutiérrez, R. Gumbs
#'
#' @export
#'
calculate_EDGE2 <- function(tree,
                            table,
                            verbose = T,
                            sort.list = FALSE){


  table <- get_extinction_prob(table, verbose = verbose)
  names(table) <- c("species", "RL.cat", "pext")

  if(isTRUE(verbose)){
    message("Calculating EDGE2 values")
  }
  N_species <- length(tree$tip.label)
  N_nodes <- tree$Nnode
  N_tot <- N_species + N_nodes

  # ensure extinction probabilities are given in same order as tree$tip.label.
  if (!identical(tree$tip.label, table$species)){
    table <- into_order(tree, table)
  }

  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }

  tree_dat <- data.frame(species = as.character(tree$tip.label),
                         TBL = NA,
                         pext = table$pext,
                         ED = NA,
                         EDGE = NA)
  tree_dat <- merge(table[,c("species", "RL.cat")], tree_dat, sort = F)



  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)

  tree <- as(tree, "phylo4")
  root <- phylobase::rootNode(tree)
  nodes <- c(root, phylobase::descendants(tree, root, "all"))

  # reorder tree components more instinctively, such that nodes are easier to find
  ord <- order(nodes)
  tree <- reorder_tree(tree, ord)
  nodes <- nodes[ord]

  tree_dat$TBL <- tree@edge.length[1:N_species]

  node_data <- data.frame(Node = 1:N_tot, Pext = rep(1, N_tot), Edge_Sum = NA)
  node_data[1:N_species, "Pext"] <- table[,"pext"]

  # assign the product of its descendant tips to each node
  for (i in c(1:length(tree@label), N_tot:(root+1))){         # for each node, beginning with tips
    anc <- tree@edge[i,1]                                   # find ancestor of node
    node_data[anc, "Pext"] <- node_data[anc, "Pext"]*node_data[i,"Pext"]   # muliply ancestor value by node "pext"
  }

  # multiply each edge.length by each pext caluclated above
  for(i in 1:length(nodes)){
    tree@edge.length[i] <- tree@edge.length[i]*node_data[i,2]
  }
  # save(tree, file ="tree.rda")

  if (is.na(tree@edge.length[root])){
    tree@edge.length[root] <- 0
  }
  node_data$Edge_Sum[root] <- tree@edge.length[root]

  # for each internal node, summate ancesteral edgelengths
  for (i in (root+1):N_tot){
    ans <- tree@edge[i,1]
    node_data$Edge_Sum[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }

  # for each tip, summate ancesteral edgelengths to find EDGE2 score
  for (i in 1:N_species){
    ans <- tree@edge[i,1]
    tree_dat$EDGE[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }

  tree_dat$ED <- tree_dat$EDGE / tree_dat$pext
  # reorder tree
  tree <- reorder_tree(tree, order(ord))

  if(isTRUE(sort.list)){
    tree_dat <- tree_dat[order(tree_dat$EDGE, decreasing = T),]
  }
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)



  edge.res <- list("list" = tree_dat, "tree" = tree, "ePDloss" = ePD.dat)

  return(edge.res)
}
