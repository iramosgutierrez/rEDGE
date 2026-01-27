
rl.data <- data.frame(
  "species" = c("Draco_drogon", "Draco_viseryi", "Draco_rheagalii"),
  "RL.cat"   = c("NT", "CR", "LC")
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


  set.seed(123)
  edge2w <- calculate_EDGE(tree, rl.data, sort.list = F, method = "EDGE2")

  expect_equal(class(edge2w), "data.frame")

  expect_message(calculate_EDGE(tree, rl.data, sort.list = F, method = "EDGE2"),
                 regexp = "Calculating EDGE2 values using Isaac extinction probabilities")

  expect_error(calculate_EDGE(tree, rl.data_broken, sort.list = F, method = "Isaac"),
               regexp = "'method' argument has to be \"EDGE1\" or \"EDGE2\"")


  rl.data_broken <- rl.data
  rl.data_broken[2,2] <- NA
  expect_error(calculate_EDGE(tree, rl.data_broken, sort.list = F, method = "EDGE2"))

})
