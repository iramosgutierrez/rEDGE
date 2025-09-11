

rl.data <- data.frame(
  "species" = c("Draco_drogon", "Draco_viseryi", "Draco_rheagalii"),
  "RL.cat"   = c("NT", "CR", "LC")
  )

tree <- ape::read.tree(text = "((Draco_rheagalii:2.5,Draco_viseryi:2.5):1,Draco_drogon:3.5);")

edge <- calculate_EDGE1(tree, rl.data)


test_that("EDGE1 runs", {

  skip_on_ci()
  skip_on_cran()

  expect(class(edge) == "data.frame", failure_message = "outupt should be a data frame")
  expect_equal(edge$ED, c(3.5, 3, 3))
  expect_equal(round(edge$EDGE, digits = 6) , c(2.197225, 1.386294, 4.158883),
         failure_message = "Unexpected EDGE1 values")

})
