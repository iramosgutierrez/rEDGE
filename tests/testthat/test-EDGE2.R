

rl.data <- data.frame(
  "species" = c("Draco_drogon", "Draco_viseryi", "Draco_rheagalii"),
  "RL.cat"   = c("NT", "CR", "LC")
)

tree <- ape::read.tree(text = "((Draco_rheagalii:2.5,Draco_viseryi:2.5):1,Draco_drogon:3.5);")

edge_df <- calculate_EDGE2(tree, rl.data, sort.list = T)
edge_ls <- calculate_EDGE2(tree, rl.data,return.all = T)

test_that("Testing EDGE2", {

expect(inherits(edge_df, "data.frame"), "Edge_df should be a data frame")
expect(inherits(edge_ls, "list"), "Edge_ls should be a list")


})






