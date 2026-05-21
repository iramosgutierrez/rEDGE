#' EDGE2 calculating function.
#'
#' @param verbose Logical. Should progress be printed or not.
#' @param return.all Logical. If TRUE, an EDGE list, tree and ePDloss list is returned. If FALSE (default), only the list is returned.
#' @param ext.prob Extinction probability associated to category, based on Isaac et al. (2007) or Mooers et al. (2008).
#'
#' @inheritParams calculate_EDGE1
#'
#' @returns A list containing 3 slots if `return.all` is TRUE. If not, only the first item (EDGE list) will be returned.
#' \itemize{
#'  \item{list: }{will output specific TBL (terminal branch length), pext (sampled probability of extinction), ED (evolutionary distinctiveness), and EDGE (EDGE index).}
#'  \item{tree: }{the input tree.}
#'  \item{ePDloss: }{complete PD (phylogenetic diversity) and ePDloss (expected PD loss). Units are as in the input tree  (generally Million years).}
#' }
#'
#' @references Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C., Baillie, J.E.M. (2007).
#' Mammals on the EDGE: Conservation Priorities Based on Threat and Phylogeny. PloS one, 2(3), e296.
#' \url{https://doi.org/10.1371/journal.pone.0000296}
#'
#' @references Mooers, A.Ø., Faith, D. P., & Maddison, W. P. (2008).
#'  Converting endangered species categories to probabilities of extinction for phylogenetic
#'  conservation prioritization. PloS one, 3(11), e3700.
#' \url{https://doi.org/10.1371/journal.pone.0003700}
#'
#' @references Gumbs, R., Gray, C. L., Böhm, M., Burfield, I. J., Couchman, O. R., Faith, D. P.,
#' Forest, F., Hoffmann, M., #' Isaac, N. J. B., Jetz, W., Mace, G. M., Mooers, A. O., Safi, K.,
#' Scott, O., Steel, M., Tucker, C. M., Pearse, W. D., Owen, N. R. & Rosindell, J. (2023).
#' The EDGE2 protocol: Advancing the prioritisation of Evolutionarily Distinct and Globally Endangered species for
#' practical conservation action. PLoS Biology, 21(2), e3001991.
#' \url{https://doi.org/10.1371/journal.pbio.3001991}
#'
#' @author I. Ramos-Gutiérrez, R. Gumbs
#'
#' @export
#'
calculate_EDGE2 <- function(tree,
                            table,
                            species.col = "species",
                            RLcat.col = "RLcat",
                            ext.prob = "Isaac",
                            verbose = TRUE,
                            sort.list  = FALSE,
                            return.all = FALSE,
                            seed = NULL){


  # check column names and rename
  if(!species.col %in% colnames(table)){
    stop("Column '", species.col, "' is not a column name in your table.\nPlease check or alternatively assign 'species' as column name in your table")
  }

  if(!RLcat.col %in% colnames(table)){
    stop("Column '", RLcat.col, "' is not a column name in your table.\nPlease check or alternatively assign 'RLcat' as column name in your table")
  }
  table <- table[,c(species.col, RLcat.col)]
  colnames(table) <-  c("species", "RLcat")

  if(is.null(seed)){
    seed <- round(runif(1, 1, 999999999))
    message(paste0("Seed has been set to: ", seed))
    }


  table <- get_extinction_prob(table, ext.prob = ext.prob, verbose = verbose, seed = seed)
  names(table) <- c("species", "RLcat", "pext")

  if(isTRUE(verbose)){
    message("Calculating EDGE2 values using ", ext.prob, " extinction probabilities")
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
  tree_dat <- merge(table[,c("species", "RLcat")], tree_dat, sort = F)



  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)

  suppressWarnings(suppressMessages(require(phylobase)))
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

  # Not tested: copied and modif. EDGE1
  # EDGEmedian <- median(tree_dat$ED)
  # tree_dat$EDGEspp <- "N"
  # tree_dat$EDGEspp[tree_dat$EDGE > EDGEmedian & tree_dat$RLcat %in% c("VU", "EN", "CR", "EW")]

  # reorder tree
  tree <- reorder_tree(tree, order(ord))

  if(isTRUE(sort.list)){
    tree_dat <- tree_dat[order(tree_dat$EDGE, decreasing = T),]
    rownames(tree_dat) <- as.character(1:nrow(tree_dat))
  }
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)



  if(isTRUE(return.all)){
  edge.res <- list("list" = tree_dat, "tree" = tree, "ePDloss" = ePD.dat)
  return(edge.res)
  }else{
  return(tree_dat)
  }
}




