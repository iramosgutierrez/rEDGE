
rl.data <- data.frame(
  "species" = c("Draco_drogon", "Draco_viseryi", "Draco_rheagalii"),
  "RLcat"   = c("NT", "CR", "LC")
)

tree <- ape::read.tree(text = "((Draco_rheagalii:2.5,Draco_viseryi:2.5):1,Draco_drogon:3.5);")


test_that("EDGE1 in wrapper", {

  skip_on_ci()
  skip_on_cran()

edge1 <- calculate_EDGE1(tree, rl.data, sort.list = F)
edgew <- calculate_EDGE(tree, rl.data, sort.list = F, method = "EDGE1")

expect_equal(edge1, edgew)

})


test_that("EDGE2 in wrapper", {

  skip_on_ci()
  skip_on_cran()



  edge2w <- calculate_EDGE(tree, rl.data, sort.list = F, method = "EDGE2", seed = 123)

  expect_equal(class(edge2w), "data.frame")

  expect_message(calculate_EDGE(tree, rl.data, sort.list = F, method = "EDGE2"),
                 regexp = "Calculating EDGE2 values using Isaac extinction probabilities")

  expect_error(calculate_EDGE(tree, rl.data_broken, sort.list = F, method = "Isaac"),
               regexp = "'method' argument has to be \"EDGE1\" or \"EDGE2\"")


  rl.data_broken <- rl.data
  rl.data_broken[2,2] <- NA
  expect_error(calculate_EDGE(tree, rl.data_broken, sort.list = F, method = "EDGE2"))

})


test_that("EDGE2 in multiphylo", {

  skip_on_ci()
  skip_on_cran()



  edge2multiphy<- calculate_EDGE(cycad.multitree[1:3], cycad.table, RLcat.col ="RL_2014", seed = 646, summarise = F)

  expect_equal(class(edge2multiphy), "list")
  expect_equal(length(edge2multiphy), 3)



})

test_that("EDGE2 several", {

  skip_on_ci()
  skip_on_cran()



  edge2multiple<- calculate_EDGE(monotreme.tree, monotreme.table,  seed = 646, summarise = T, n.iter = 10)

  expect_equal(inherits(edge2multiple, "data.frame"), TRUE)
  expect_equal(dim(edge2multiple), c(length(monotreme.tree$tip.label), 13))

})

test_that("EDGE2 several and mutiphylo", {

  skip_on_ci()
  skip_on_cran()


  edgeall<- calculate_EDGE(cycad.multitree[1:3], cycad.table, RLcat.col ="RL_2014", seed = 646, summarise = F, n.iter = 4)
  expect_equal(inherits(edgeall, "list"), TRUE)
  expect_equal(nrow(dplyr::bind_rows(edgeall)), 12*length(cycad.multitree[[1]]$tip.label))
})
