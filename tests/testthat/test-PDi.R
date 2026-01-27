
rl.data <- rEDGE::crocodile.table
trees <- rEDGE::crocodile.trees
tree <- trees[1]


test_that("PDindicator", {

  skip_on_ci()
  skip_on_cran()

  pdi <- calculate_PD_indicator(trees[1:3],
                                rl.data,
                                time.cols = c("RL.2023", "RL.2024"),

                                verbose = FALSE, parallelize = F,
                                seed = 123)

  expect_equal(class(pdi), "list")
  expect_equal(names(pdi), c("RL.2023", "RL.2024"))
})


test_that("different extprob", {

  skip_on_ci()
  skip_on_cran()

  expect_message(calculate_PD_indicator(trees[1:3], rl.data,
                                        time.cols = c("RL.2023", "RL.2024"),

                                verbose = T,
                                seed = 123),

                 regexp = "Calculating EDGE2 values using Isaac extinction probabilities"
  )

  expect_message(calculate_PD_indicator(trees[1:3], rl.data,
                                        time.cols = c("RL.2023", "RL.2024"),

                                        verbose = TRUE, seed = 123, ext.prob = "IUCN50"),

                 regexp = "Calculating EDGE2 values using IUCN50 extinction probabilities"
  )


})

