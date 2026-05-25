
rl.data <- rEDGE::cycad.table
tree <- rEDGE::cycad.tree
trees <- rEDGE::cycad.multitree



test_that("PDindicator", {

  skip_on_ci()
  skip_on_cran()

  pdi <- calculate_PD_indicator(tree,
                                rl.data,
                                time.cols = c("RL_2003", "RL_2014"),

                                verbose = FALSE, parallelize = F,
                                seed = 123)

  expect_equal(class(pdi), "list")
  expect_equal(names(pdi), time.cols)
})


test_that("different extprob", {

  skip_on_ci()
  skip_on_cran()

  expect_message(calculate_PD_indicator(tree, cycad.table,
                                        time.cols = c("RL_2003", "RL_2014"),

                                verbose = T,
                                seed = 123),

                 regexp = "Calculating EDGE2 values using Isaac extinction probabilities"
  )

  expect_message(calculate_PD_indicator(trees[1:3], cycad.table,
                                        time.cols = c("RL_2003", "RL_2014"),

                                        verbose = TRUE, seed = 123, ext.prob = "IUCN50"),

                 regexp = "Calculating EDGE2 values using IUCN50 extinction probabilities"
  )


})

