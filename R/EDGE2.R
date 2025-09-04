

calculate_EDGE2 <- function(tree, table, verbose = F){

  # require(phylobase)
  # require(data.table)

  table <- get_extinction_prob(table)
  names(table) <- c("species", "RL.cat", "pext")

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

  tree_dat <- data.frame(Species = as.character(tree$tip.label),
                         TBL = NA,
                         pext = table$pext,
                         ED = NA,
                         EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),
                        ePDloss = NA)


  # tree <- suppressWarnings(as(tree, "phylo4"))
  # root <- rootNode(tree)
  # nodes <- c(root, descendants(tree, root, "all"))
  #### root <- ape::Ntip(tree)+1


  # # reorder tree components more instinctively, such that nodes are easier to find
  # ord <- order(nodes)
  # tree <- reorder_tree(tree, ord)
  # nodes <- nodes[ord]

  tree_table <- cbind(tree$edge, tree$edge.length)
  colnames(tree_table) <- c("parent", "node", "length")

  tipnames <- data.frame("node" = 1:N_species,
                         "species" = tree$tip.label)
  nodenames <- data.frame("node" = N_species+1:N_tot,
                         "species" = NA)

  tree_table <- merge(tree_table, rbind(tipnames, nodenames) )

  tree_dat$TBL <- tree@edge.length[1:N_species]

  node_data <- data.frame(Node = 1:N_tot, Pext = rep(1, N_tot), Edge_Sum = NA)
  node_data[1:N_species, "Pext"] <- table[,"pext"]

  # assign the product of its descendant tips to each node
  for (i in c(1:length(tree@label), N_tot:(root+1))){         # for each node, beginning with tips
    anc <- tree@edge[i,1]                                     # find ancestor of node
    node_data[anc, 2] <- node_data[anc, 2]*node_data[i,2]     # muliply ancestor value by node "pext"
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

  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}
