
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



cat_pext <- function(){
  cat.pext <- data.frame(rl.cat = rev(c("CR","EN","VU","NT","LC")) ,
                         pext = rev(c(0.97, 0.97/2, 0.97/4,0.97/8,0.97/16)))
  return(cat.pext)
}


# GE2 calculation #
# generate 1,000,000 GE2 values distributed across all Red List categories
# and sample to a distribution of GE2 values to each RL category to capture uncertainty
# this distribution of GE2 is then used to derive the GE2 for each species for each iteration of EDGE2 calculation
# a random value of GE2, corresponding to the appropriate RL category, can be assigned from the output to each species

create_pext_by_cat <- function(n = 1000000){

  iucn <- sample(1:5, size=n, replace=TRUE) # 1:5 forces to be just 5 categories!
  data <- data.frame(species=c(1:n), pext=cat_pext()$pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  pext <- cat_pext()
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext.tap <- c(0, pext$pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4; rank.quu <- rank^5
  model <- lm(pext.tap ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2; data$rank.cub <- data$rank^3; data$rank.qu <- data$rank^4; data$rank.quu <- data$rank^5
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  pext.LC <- data.frame(RL.cat = "LC", pext =data$rank.pext[data$pext == pext.tap[2]])
  pext.NT <- data.frame(RL.cat = "NT", pext =data$rank.pext[data$pext == pext.tap[3]])
  pext.VU <- data.frame(RL.cat = "VU", pext =data$rank.pext[data$pext == pext.tap[4]])
  pext.EN <- data.frame(RL.cat = "EN", pext =data$rank.pext[data$pext == pext.tap[5]])
  pext.CR <- data.frame(RL.cat = "CR", pext =data$rank.pext[data$pext == pext.tap[6]])
  pext.obj <- rbind(pext.CR,pext.EN, pext.VU, pext.NT, pext.LC)

  pext.obj.sample <- NULL



  #    I don't get this!
  #
  #

  # for(i in unique(pext.obj$RL.cat)){
  #
  #   a <- pext.obj$pext[pext.obj$RL.cat == i]
  #   if(median(a) < pext$pext[pext$rl.cat == i]){
  #     while(median(a) < pext$pext[pext$rl.cat == i]){
  #       a <- a[-sample(c(1:length(a)),50,prob = rev(a))]
  #     }
  #   }else{
  #     while(median(a) > pext$pext[pext$rl.cat == i]){
  #       a <- a[-sample(c(1:length(a)),50,prob = a)]
  #     }
  #   }
  #   pext.obj.sample <- rbind(pext.obj.sample,data.frame(RL.cat = i, pext = a))
  #   progressbar(which(unique(pext.obj$RL.cat) == i), length(unique(pext.obj$RL.cat)),
  #               msg = "Sampling extinction probabilities...")
  # }
  # return(pext.obj.sample)
  return(pext.obj)
}



#' Probability of Extinction
#'
#'
#'
get_extinction_prob <- function(table, verbose = T){

  if(!inherits(table, c("data.frame", "tibble"))){
    stop("table object should be a tibble or data.frame")
  }

  if(!all(colnames(table) == c("species", "RL.cat"))){
    stop("Column names of 'table' should be 'species' and RL.cat'")
  }

  cat_pext_table <- create_pext_by_cat()


  for(sp in table$species){

    if(isTRUE(verbose)){
      progressbar(which(table$species == sp),
                  length(table$species),
                  msg = "Calculating extinction probabilites...")
    }
    cat.i <- table$RL.cat[table$species==sp]
    if(cat.i =="NE"){pext.i <- runif(1,0.0001, 0.9999)}else
      if(cat.i =="EX"){
        cat.i =="CR"
      }else{
        pext.i <- sample(cat_pext_table$pext[cat_pext_table$RL.cat==cat.i], size = 1)
      }
    table$pext[table$species== sp] <- pext.i
  }
  return(table)
}


