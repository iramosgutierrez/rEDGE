
progressbar <- function (curr.iter, tot.iter, ini.iter = 1, units = "mins",
                         msg = NULL){
  curr.iter <- curr.iter - ini.iter + 1
  tot.iter <- tot.iter - ini.iter + 1
  if (units == "secs") {
    d <- 0
  }
  else if (units == "hours") {
    d <- 2
  }
  else {
    d <- 1
  }
  if (curr.iter == 1 & !is.null(msg)) {
    cat(msg, "\n")
  }
  if (curr.iter == 1) {
    st <<- Sys.time()
    cat(paste0("0%       25%       50%       75%       100%",
               "\n", "|---------|---------|---------|---------|",
               "\n"))
  }
  v <- seq(from = 0, to = 40, by = 40/tot.iter)
  v <- diff(ceiling(v))
  v <- cumsum(v)
  txt <- strrep("*", times = v[curr.iter])
  txt <- stringr::str_pad(txt, width = 45, side = "right",
                          pad = " ")
  ct <- Sys.time()
  et <- as.numeric(difftime(ct, st, units = units))/curr.iter *
    (tot.iter - curr.iter)
  et <- round(et, digits = d)
  txt.end <- paste0(txt, "ETC: ", et, " ", units)
  if (curr.iter == ini.iter) {
    txt.end <- paste0(txt, "ETC: ")
    maxnchar <<- nchar(txt.end)
  }
  if (curr.iter == tot.iter) {
    txt.end <- paste0("*", txt, "DONE")
  }
  if (nchar(txt.end) > maxnchar) {
    maxnchar <<- nchar(txt.end)
  }
  txt.end <- stringr::str_pad(txt.end, width = maxnchar, side = "right",
                              pad = " ")
  cat("\r")
  cat(txt.end)
  if (curr.iter == tot.iter) {
    cat("\n")
  }
  if (curr.iter == tot.iter) {
    rm(list = c("st", "maxnchar"), envir = .GlobalEnv)
  }
}



cat_pext <- function(ext.prob = "Isaac"){
  if(!ext.prob %in% c("Isaac", "IUCN50", "IUCN100", "IUCN500")){
    stop("ext.prob should be one of: 'Isaac' / 'IUCN50' / 'IUCN100' / 'IUCN500'")
    }

  if(ext.prob == "Isaac"){
  cat.pext <- data.frame(rl.cat = rev(c("CR", "EN",   "VU",   "NT",   "LC")) ,
                         pext   = rev(c(0.97, 0.97/2, 0.97/4, 0.97/8, 0.97/16)))
  }

  if(ext.prob == "IUCN50"){
    cat.pext <- data.frame(rl.cat = rev(c("CR", "EN",  "VU",   "NT",   "LC")) ,
                           pext   = rev(c(0.97,  0.42,  0.05,   0.004,  0.00005)))
  }

  if(ext.prob == "IUCN100"){
    cat.pext <- data.frame(rl.cat = rev(c("CR",  "EN",   "VU",   "NT",   "LC")) ,
                           pext   = rev(c( 0.999, 0.667,  0.1,    0.01,   0.0001)))
  }

  if(ext.prob == "IUCN500"){
    cat.pext <- data.frame(rl.cat = rev(c("CR", "EN",   "VU",   "NT",   "LC")) ,
                           pext   = rev(c( 1,    0.996,  0.39,   0.02,   0.0005)))
  }

  return(cat.pext)
}


# GE2 calculation #
# generate 1,000,000 GE2 values distributed across all Red List categories
# and sample to a distribution of GE2 values to each RL category to capture uncertainty
# this distribution of GE2 is then used to derive the GE2 for each species for each iteration of EDGE2 calculation
# a random value of GE2, corresponding to the appropriate RL category, can be assigned from the output to each species


make_pext_curve <- function(ext.prob){

  x_pts <- 1:5
  y_pts <- cat_pext(ext.prob)$pext

  # avoid logit blow-up at 0/1
  eps <- 1e-12
  y_pts <- pmin(pmax(y_pts, eps), 1 - eps)

  logit_y <- log(y_pts / (1 - y_pts))

  phi <- splinefun(x_pts, logit_y, method = "monoH.FC")

  function(x){
    plogis(phi(x))  # back to probability space
  }
}

create_pext_by_cat <- function(n = 1000000, ext.prob){


  pext.dist <- data.frame(RL.num=seq(0, 6, length.out = n)) |>
    dplyr::mutate(pext = make_pext_curve(ext.prob)(RL.num))

    pext.LC <- data.frame(RLcat = "LC", pext = pext.dist$pext[pext.dist$RL.num >= 0.0 & pext.dist$RL.num < 1.5])
    pext.NT <- data.frame(RLcat = "NT", pext = pext.dist$pext[pext.dist$RL.num >= 1.5 & pext.dist$RL.num < 2.5])
    pext.CD <- data.frame(RLcat = "CD", pext = pext.dist$pext[pext.dist$RL.num >= 1.5 & pext.dist$RL.num < 2.5]) # Conservation dependent ~ NT
    pext.VU <- data.frame(RLcat = "VU", pext = pext.dist$pext[pext.dist$RL.num >= 2.5 & pext.dist$RL.num < 3.5])
    pext.EN <- data.frame(RLcat = "EN", pext = pext.dist$pext[pext.dist$RL.num >= 3.5 & pext.dist$RL.num < 4.5])
    pext.CR <- data.frame(RLcat = "CR", pext = pext.dist$pext[pext.dist$RL.num >= 4.5 & pext.dist$RL.num < 5.5])
    pext.EW <- data.frame(RLcat = "EW", pext = pext.dist$pext[pext.dist$RL.num >= 5.5 & pext.dist$RL.num < 6.0]) #EW ~ CR

    pext.obj <- rbind(pext.LC, pext.NT, pext.CD, pext.VU, pext.EN, pext.CR, pext.EW)
    return(pext.obj)
}

create_pext_by_cat_LM_defunct <- function(n = 1000000, ext.prob){

  iucn <- sample(1:5, size=n, replace=TRUE) # 1:5 forces to be just 5 categories!
  data <- data.frame(species=c(1:n), pext=cat_pext(ext.prob)$pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  pext <- cat_pext(ext.prob)
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext.tap <- c(0, pext$pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4; rank.quu <- rank^5
  model <- lm(pext.tap ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2; data$rank.cub <- data$rank^3; data$rank.qu <- data$rank^4; data$rank.quu <- data$rank^5
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  pext.LC <- data.frame(RLcat = "LC", pext =data$rank.pext[data$pext == pext.tap[2]])
  pext.NT <- data.frame(RLcat = "NT", pext =data$rank.pext[data$pext == pext.tap[3]])
  pext.CD <- data.frame(RLcat = "CD", pext =data$rank.pext[data$pext == pext.tap[3]]) # Conservation dependent ~ NT
  pext.VU <- data.frame(RLcat = "VU", pext =data$rank.pext[data$pext == pext.tap[4]])
  pext.EN <- data.frame(RLcat = "EN", pext =data$rank.pext[data$pext == pext.tap[5]])
  pext.CR <- data.frame(RLcat = "CR", pext =data$rank.pext[data$pext == pext.tap[6]])
  pext.EW <- data.frame(RLcat = "EW", pext =data$rank.pext[data$pext == pext.tap[6]]) #EW ~ CR
  pext.obj <- rbind(pext.CR,pext.EN, pext.VU, pext.NT, pext.LC)

  pext.obj.sample <- NULL



  #    I don't get this!
  #
  #

  # for(i in unique(pext.obj$RLcat)){
  #
  #   a <- pext.obj$pext[pext.obj$RLcat == i]
  #   if(median(a) < pext$pext[pext$RLcat == i]){
  #     while(median(a) < pext$pext[pext$RLcat == i]){
  #       a <- a[-sample(c(1:length(a)),50,prob = rev(a))]
  #     }
  #   }else{
  #     while(median(a) > pext$pext[pext$rl.cat == i]){
  #       a <- a[-sample(c(1:length(a)),50,prob = a)]
  #     }
  #   }
  #   pext.obj.sample <- rbind(pext.obj.sample,data.frame(RLcat = i, pext = a))
  #   progressbar(which(unique(pext.obj$RLcat) == i), length(unique(pext.obj$RLcat)),
  #               msg = "Sampling extinction probabilities...")
  # }
  # return(pext.obj.sample)
  return(pext.obj)
}




get_extinction_prob <- function(table, ext.prob = ext.prob, verbose = T, seed = seed){

  if(!inherits(table, c("data.frame", "tibble"))){
    stop("table object should be a tibble or data.frame")
  }

  if(!all(colnames(table) == c("species", "RLcat"))){
    stop("Column names of 'table' should be 'species' and 'RLcat'")
  }

  set.seed(seed*2)
  cat_pext_table <- create_pext_by_cat(ext.prob = ext.prob)

  table$pext <- NA
  for(sp in table$species){

    if(isTRUE(verbose)){
      progressbar(which(table$species == sp),
                  length(table$species),
                  msg = "Calculating extinction probabilites...")
    }
    cat.i <- table$RLcat[table$species==sp]
    if(cat.i %in% c("NE", "DD")){
      # pext.i <- runif(1,0.0001, 0.9999)  # this makes DD and NE get ~ 0.5 pext in average
      set.seed(seed + which(table$species == sp))
      pext.i <- sample(cat_pext_table$pext[cat_pext_table$pext < 0.999], size = 1)
    }else{
      if(cat.i %in% c("EX", "EW")){cat.i <- "CR"}
      if(cat.i =="CD"){cat.i <-"NT"}

      set.seed(seed + 2*which(table$species == sp))
      pext.i <- sample(cat_pext_table$pext[cat_pext_table$RLcat==cat.i], size = 1)
      }
    table$pext[table$species== sp] <- pext.i
  }
  return(table)
}


IQR <- function(x, na.rm = T){ #Inter-quartilic range
  return(stats::quantile(x, 0.75, na.rm = na.rm) - stats::quantile(x, 0.25, na.rm = na.rm))
}


# For EDGE2 function
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

# for plotting
iucn_palette <- c("NA" = "#C1B5A5",
                  "NE" = "#FFFFFF",
                  "DD" = "#D1D1C6",
                  "LC" = "#60C659",
                  "NT" = "#CCE226",
                  "VU" = "#F9E814",
                  "EN" = "#FC7F3F",
                  "CR" = "#D81E05",
                  "RE" = "#9B4F96",
                  "EW" = "#542344",
                  "EX" = "#000000")
# levels(iucn_palette) <- c("NA",  "NE", "DD", "LC", "NT",  "VU", "EN", "CR", "RE",  "EW", "EX")
